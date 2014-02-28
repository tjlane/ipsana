#!/usr/bin/env python

import os
import sys
import time
import argparse

import numpy as np
import psana 

import zmq

from psdata import ImageData
from psdata import XYPlotData

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
    '-g',
    '--geometry',
    type=str,
    default='default_geometry.npy',
    help='A numpy array on disk specifying the q-values for each bin'
)

parser.add_argument(
    '-m',
    '--mask',
    type=str,
    default='default_mask.npy',
    help='A numpy array on disk specifying the masked pixels'
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
print "Loading geometry from:        %s" % args.geometry
geometry = np.abs(np.load(args.geometry))
print "Loading pixel mask from:      %s" % args.mask
mask = np.load(args.mask)

# Define experiment, run. Shared memory syntax (for SXR): shmem=0_1_SXR.0:stop=no
if args.run >= 0:
    print "\nReading run #: %d" % args.run
    ds = psana.DataSource('exp=xppb0114:run=%d' % args.run)
elif args.run == -1:
    print "\nLocking on to shared memory..."
    ds = psana.DataSource('shmem=0_1_XPP.0:stop=no')
else:
    raise ValueError('`run` parameter must either be an int corresponding to a '
                     'run, or the keyword "online"')

cspad_src  = psana.Source('DetInfo(XppGon.0:Cspad.0)')
evr_src    = psana.Source('DetInfo(NoDetector.0:Evr.0)')

plot_update_frequency = 10


# ------ END CONFIG ---------

epics = ds.env().epicsStore()
calib = ds.env().calibStore()

laser_on  = np.zeros((32, 185, 194*2))
laser_off = np.zeros((32, 185, 194*2))

shot_index = 0

n_laser_on_shots  = 0
n_laser_off_shots = 0

shot_result = 'no data'

for i,evt in enumerate(ds.events()):
    
    EVR_code = evt.get(psana.EvrData.DataV3, evr_src)
    if EVR_code is None:
        continue
        
    cspad = evt.get(psana.CsPad.DataV2, cspad_src)
    #corrected_cspad = evt.get(psana.CsPad.DataV2, cspad_src, 'calibrated_ndarr')

    # extract useful data
    raw_image       = np.vstack([ cspad.quads(i).data() for i in range(4) ])
    #corrected_image = np.vstack([ corrected_cspad.quads(i).data() for i in range(4) ])
    corrected_image = raw_image.copy()

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


    # === compute radial averages ===
    bins, evt_rad     = radial_average(corrected_image, geometry, mask)
    bins, avg_rad_on  = radial_average(laser_on, geometry, mask)
    bins, avg_rad_off = radial_average(laser_off, geometry, mask)

    # ==============================
    # PUBLISH DATA VIA ZMQ    
    # ==============================
    rad_avg_plot = XYPlotData(n_laser_off_shots + n_laser_on_shots,
                               "I(q) Laser On & Off",
                              bins, evt_rad)
    socket.send("radavg", zmq.SNDMORE)
    socket.send_pyobj(rad_avg_plot)
        
    print "Run: %d | Shot %d | %s" % (args.run, shot_index, shot_result)
    shot_index += 1

