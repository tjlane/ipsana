#!/usr/bin/env python

import os
import sys
import time
import argparse

import numpy as np
import psana 

from utilities import *

# ARGUMENT PARSING
parser = argparse.ArgumentParser(description='Psana plot server application')

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
print '\t\t MINITTI BEAMTIME TIMETOOL ANALYSIS -- PEW PEW PEW'
print ''
print '*' * 80
print ''


# PSANA CONFIG FILE - this must be set before a datasource is created
basedir = os.path.split(os.path.abspath( __file__ ))[0]
config_fn = os.path.join(basedir, "minitti.cfg")
psana.setConfigFile(config_fn)
psana.setOption('psana.l3t-accept-only',0)
print "Loading psana config file:    %s" % config_fn

# Aquire the geometry and mask
geometry_filename = '/reg/neh/home2/tjlane/analysis/xppb0114/geometries/v1/v2_q_geom.npy'
print "Loading geometry from:        %s" % geometry_filename
geometry = np.load(geometry_filename)

mask_filename = '/reg/neh/home2/tjlane/analysis/xppb0114/geometries/v1/mask.npy'
print "Loading pixel mask from:      %s" % mask_filename
mask = np.load(mask_filename)

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

n_q_bins = 101
ra = RadialAverager(geometry, mask, n_bins=n_q_bins)
tt_hist = TTHistogram(-0.4, 0.4, 101, n_q_bins)

for i,evt in enumerate(ds.events()):
    
    EVR_code = evt.get(psana.EvrData.DataV3, evr_src)
    if EVR_code is None:
        continue
        
    corrected_image = evt.get(psana.ndarray_float32_3, cspad_src, 'calibrated_ndarr')
    assert corrected_image.shape == (32, 185, 194*2)

    xfel_status, uv_status = pumpprobe_status(EVR_code)
    
    # === filter for bad shots ===
    if corrected_image.sum() < args.threshold:
        print "-- bad shot :: intensity too low"
        continue

    # === extract timetool stamp ===
    tau = epics.value('TTSPEC:FLTPOS_PS')

    # === accumulate the results ===
    if (xfel_status and uv_status):       # good data, w/ pump
        shot_result = 'uv & xfel on'
        n_laser_on_shots += 1
        update_average(n_laser_on_shots, laser_on, corrected_image)
        
    else:                                 # no scattering data...
        shot_result = 'uv off'

    # === ANALYZE ===
    q_values, evt_rad = ra(corrected_image)
    n_evt_rad = normalize(q_values, evt_rad)
    tt_hist.add_data(evt_rad, tau)
        
    print "Run: %d | Shot %d | %s" % (args.run, shot_index, shot_result)
    shot_index += 1

print tt_hist._histogram.sum()
tt_hist.save('/reg/neh/home2/tjlane/analysis/xppb0114/tt_hists/r%d_tthist.h5' % args.run)
