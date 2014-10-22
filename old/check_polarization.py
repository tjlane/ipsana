#!/usr/bin/env python

import os
import sys
import time
import h5py

import numpy as np
import psana 

from psdata import ImageData
from psdata import IqPlotData

from utilities import *

# =============================================
s_runs = [182, 183]
p_runs = [193, 194]
# ==============================================

# PSANA CONFIG FILE - this must be set before a datasource is created
basedir = os.path.split(os.path.abspath( __file__ ))[0]
config_fn = os.path.join(basedir, "minitti.cfg")
psana.setConfigFile(config_fn)
psana.setOption('psana.l3t-accept-only',0)
print "Loading psana config file:    %s" % config_fn

# Aquire the geometry and mask
geometry_filename = '/reg/neh/home2/tjlane/analysis/xppb0114/geometries/v2/q_geom.npy'
print "Loading geometry from:        %s" % geometry_filename
geometry = np.load(geometry_filename).reshape(32,185,388)

mask_filename = '/reg/neh/home2/tjlane/analysis/xppb0114/geometries/v2/mask_v2.npy'
print "Loading pixel mask from:      %s" % mask_filename
mask = np.load(mask_filename).reshape(32,185,388)


cspad_src  = psana.Source('DetInfo(XppGon.0:Cspad.0)')
evr_src    = psana.Source('DetInfo(NoDetector.0:Evr.0)')

update_frequency = 100


# ------ END CONFIG ---------

p_polarized = np.zeros((32, 185, 194*2))
s_polarized = np.zeros((32, 185, 194*2))

n_q_bins = 100
ra = RadialAverager(geometry, mask, n_bins=n_q_bins)


def sum_images(run):

  ds = psana.DataSource('exp=xppb0114:run=%d' % run)
  epics = ds.env().epicsStore()
  calib = ds.env().calibStore()

  shot_index = 0
  good_shot_index = 0
  average_image = np.zeros((32, 185, 194*2))

  for i,evt in enumerate(ds.events()):
    
    EVR_code = evt.get(psana.EvrData.DataV3, evr_src)
    if EVR_code is None:
        continue
        
    cspad = evt.get(psana.CsPad.DataV2, cspad_src)
    corrected_image = evt.get(psana.ndarray_float32_3, cspad_src, 'calibrated_ndarr')

    assert corrected_image.shape == (32, 185, 194*2)

    xfel_status, uv_status = pumpprobe_status(EVR_code)
    
    # === filter for bad shots ===
    
    # filter 
    if corrected_image.sum() < 0.0:
        print "-- bad shot :: intensity too low"
        continue

    # === extract timetool stamp ===

    tau = epics.value('TTSPEC:FLTPOS_PS')
    angle_shift = epics.value('LAS:FS3:REG:Angle:Shift:rd')

    # === accumulate the results ===

    if (xfel_status and uv_status):       # good data, w/ pump
        shot_result = 'uv & xfel on'
        good_shot_index += 1
        update_average(shot_index, average_image, corrected_image)
        
    elif (xfel_status and not uv_status): # good data, no pump
        shot_result = 'uv off, xfel on'
        
    else:                                 # no scattering data...
        shot_result = 'xfel off'

    print "Run: %d | Shot %d | %s" % (run, shot_index, shot_result)

    shot_index += 1

  return average_image


for i,run in enumerate(p_runs):
    update_average(i, p_polarized, sum_images(run))

for i,run in enumerate(s_runs):
    update_average(i, s_polarized, sum_images(run))
    
print '\nSaving final results...'

hdf = h5py.File('/reg/neh/home2/tjlane/analysis/xppb0114/polarization/cspad_sp_difference_image.h5')
hdf['/p_polarized'] = p_polarized
hdf['/s_polarized'] = s_polarized
hdf.close()

