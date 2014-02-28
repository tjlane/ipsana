#!/usr/bin/env python

import os
import sys
import time
import argparse

import numpy as np
import psana 

import zmq

from psdata import ImageData
from psdata import IqPlotData

from utilities import *

# ARGUMENT PARSING
parser = argparse.ArgumentParser(description='Psana plot server application')

parser.add_argument(
    '-p',
    '--port',
    metavar='PORT',
    type=int,
    default=5556,
    help='the tcp port to use on the server'
)

parser.add_argument(
    '-r',
    '--run',
    type=int,
    default=-1,
    help='Which run to analyze, -1 for live stream',
)

parser.add_argument(
    '-t',
    '--threshold',
    type=int,
    default=0,
    help='An intensity threshold value',
)


args = parser.parse_args()

print ''
print '*' * 80
print ''
print '\t\t MINITTI BEAMTIME SERVER -- PEW PEW PEW'
print ''
print '*' * 80
print ''

# ZMQ SETUP 
context = zmq.Context()
socket = context.socket(zmq.PUB)
socket.setsockopt(zmq.SNDHWM, 10)
socket.bind("tcp://*:%d" % args.port)
print "Broadcasting via ZMQ on port: %d" % args.port

# PSANA CONFIG FILE - this must be set before a datasource is created
basedir = os.path.split(os.path.abspath( __file__ ))[0]
config_fn = os.path.join(basedir, "minitti.cfg")
psana.setConfigFile(config_fn)
psana.setOption('psana.l3t-accept-only',0)
print "Loading psana config file:    %s" % config_fn

# Aquire the geometry and mask
geometry_filename = '/reg/neh/home2/tjlane/analysis/xppb0114/geometries/v1/v2_q_geom.npy'
print "Loading geometry from:        %s" % geometry_filename
geometry = np.load(geometry_filename).reshape(32,185,388)

mask_filename = '/reg/neh/home2/tjlane/analysis/xppb0114/geometries/v1/mask.npy'
print "Loading pixel mask from:      %s" % mask_filename
mask = np.load(mask_filename).reshape(32,185,388)

# Define experiment, run. Shared memory syntax (for SXR): shmem=0_1_SXR.0:stop=no
if args.run >= 0:
    print "\nReading run #: %d" % args.run
    ds = psana.DataSource('exp=xppb0114:run=%d' % args.run)
elif args.run == -1:
    print "\nLocking on to shared memory..."
    ds = psana.DataSource('shmem=1_1_psana_XPP.0:stop=no')
else:
    raise ValueError('`run` parameter must either be an int corresponding to a '
                     'run, or the keyword "online"')

cspad_src  = psana.Source('DetInfo(XppGon.0:Cspad.0)')
evr_src    = psana.Source('DetInfo(NoDetector.0:Evr.0)')

update_frequency = 100


# ------ END CONFIG ---------

epics = ds.env().epicsStore()
calib = ds.env().calibStore()

laser_on  = np.zeros((32, 185, 194*2))
laser_off = np.zeros((32, 185, 194*2))

shot_index = 0

n_laser_on_shots  = 0
n_laser_off_shots = 0

shot_result = 'no data'
t0 = time.time()
s0 = 0 # shot counter

ra = RadialAverager(geometry, mask)

for i,evt in enumerate(ds.events()):
    
    EVR_code = evt.get(psana.EvrData.DataV3, evr_src)
    if EVR_code is None:
        continue
        
    cspad = evt.get(psana.CsPad.DataV2, cspad_src)
    corrected_image = evt.get(psana.ndarray_float32_3, cspad_src, 'calibrated_ndarr')
    raw_image       = np.vstack([ cspad.quads(i).data() for i in range(4) ])

    assert raw_image.shape == (32, 185, 194*2)
    assert corrected_image.shape == (32, 185, 194*2)

    xfel_status, uv_status = pumpprobe_status(EVR_code)
    
    # === filter for bad shots ===
    
    # filter 
    if corrected_image.sum() < args.threshold:
        print "-- bad shot :: intensity too low"
        continue
    
    # === accumulate the results ===

    if (xfel_status and uv_status):       # good data, w/ pump
        shot_result = 'uv & xfel on'
        n_laser_on_shots += 1
        update_average(n_laser_on_shots, laser_on, corrected_image)
        
    elif (xfel_status and not uv_status): # good data, no pump
        shot_result = 'uv off, xfel on'
        n_laser_off_shots += 1
        update_average(n_laser_off_shots, laser_off, corrected_image)
        
    else:                                 # no scattering data...
        shot_result = 'xfel off'

    if ((n_laser_on_shots + n_laser_off_shots) % update_frequency == 0) \
        and (n_laser_on_shots > 0):

        # === compute radial averages ===
        bins1, evt_rad     = ra(corrected_image)
        bins2, avg_rad_on  = ra(laser_on)
        bins3, avg_rad_off = ra(laser_off)

        # normalize from 0.5 to 3.5 inv ang
        n_evt_rad     = normalize(bins1, evt_rad)
        n_avg_rad_on  = normalize(bins2, avg_rad_on)
        n_avg_rad_off = normalize(bins3, avg_rad_off)

        # --- any data massaging for laser on/off ---
        laser_onoff_diff = n_avg_rad_on - n_avg_rad_off
        laser_onoff_percdiff = laser_onoff_diff / n_avg_rad_on

        # ==============================
        # PUBLISH DATA VIA ZMQ    
        # ==============================
        rad_evt_plot = IqPlotData(n_laser_off_shots + n_laser_on_shots,
                                "I(q) Instantaneous",
                                 bins1, evt_rad)
        socket.send("radevt", zmq.SNDMORE)
        socket.send_pyobj(rad_evt_plot)


        plotdata = np.hstack((avg_rad_on[:,None], avg_rad_off[:,None]))
        rad_avg_plot = IqPlotData(n_laser_off_shots + n_laser_on_shots,
                                "I(q) Laser On & Off Averages",
                                 bins2, plotdata)
        socket.send("radavg", zmq.SNDMORE)
        socket.send_pyobj(rad_avg_plot)


        rad_diff_plot = IqPlotData(n_laser_off_shots + n_laser_on_shots,
                                "I(q) Laser On/Off Difference",
                                 bins3, laser_onoff_diff)
        socket.send("raddiff", zmq.SNDMORE)
        socket.send_pyobj(rad_diff_plot)


        rad_diffperc_plot = IqPlotData(n_laser_off_shots + n_laser_on_shots,
                                "I(q) Laser On/Off Difference (%)",
                                 bins3, laser_onoff_percdiff)
        socket.send("radperc", zmq.SNDMORE)
        socket.send_pyobj(rad_diffperc_plot)


        delta_t = time.time() - t0
        t0 = time.time()

        shot_rate = float(shot_index - s0) / delta_t
        s0 = shot_index

        print '\t--> sending results (%d Hz)' % (shot_rate)
        
    print "Run: %d | Shot %d | %s" % (args.run, shot_index, shot_result)
    shot_index += 1


print '\nSaving final results...'

np.save('/reg/neh/home2/tjlane/analysis/xppb0114/averages/r%d_laser_on_avg.npy' % args.run, avg_rad_on)
np.save('/reg/neh/home2/tjlane/analysis/xppb0114/averages/r%d_laser_off_avg.npy' % args.run, avg_rad_off)


