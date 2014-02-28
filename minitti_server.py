#!/usr/bin/env python

import os
import sys
import time

import numpy as np
import psana 

import zmq

from psdata import ImageData
from psdata import XYPlotData

import utilities

################################################################################
# ZMQ SETUP 
context = zmq.Context()
socket = context.socket(zmq.PUB)
socket.setsockopt(zmq.SNDHWM, 10)
socket.bind("tcp://*:%d" % 12322)
print "Broadcasting ZMQ on port 12322"
################################################################################

# ------- CONFIG -------
# PSANA CONFIG FILE - this must be set before a datasource is created
basedir = os.path.split(os.path.abspath( __file__ ))[0]
config_fn = os.path.join(basedir, "minitti.cfg")
psana.setConfigFile(config_fn)
psana.setOption('psana.l3t-accept-only',0) 

# Define experiment, run. Shared memory syntax (for SXR): shmem=0_1_SXR.0:stop=no
ds = psana.DataSource('exp=xppi0613:run=275') 

cspad_src  = psana.Source('DetInfo(XppGon.0:Cspad.0)')
evr_src    = psana.Source('DetInfo(NoDetector.0:Evr.0)')
FEEGas_src = psana.Source('BldInfo(FEEGasDetEnergy)')

status_frequency      = 100         # print a status every on the Nth shot
plot_update_frequency = 100         # carry out radial integration and update plots every N good shots


# ------ END CONFIG ---------

epics = ds.env().epicsStore()
calib = ds.env().calibStore()

for i,evt in enumerate(ds.events()):
    
    EVR_code = evt.get(psana.EvrData.DataV3, evr_src)
    if EVR_code is None:
        continue
        
    cspad = evt.get(psana.CsPad.DataV2, cspad_src)
    corrected_cspad = evt.get(psana.CsPad.DataV2, src, 'calibrated_ndarr')

    # extract useful data
    raw_image       = np.vstack([ cspad.quads(i).data() for i in range(4) ])
    corrected_image = np.vstack([ corrected_cspad.quads(i).data() for i in range(4) ])
    assembled_image = evt.get(psana.ndarray_float32_2, cspad_src, 'assembled_image')

    print raw_image.shape
    print corrected_iamge.shape
    print assembled_image.shape

    xfel_status, uv_status = pumpprob_status(EVR_code)



    # ==============================
    # PUBLISH DATA VIA ZMQ    
    # ==============================
    avgImage = ImageData(goodEvents, "AVERAGE OF ALL CALIBRATED CSPAD", avgCspadAll)
    socket.send("avgcspad", zmq.SNDMORE)
    socket.send_pyobj(avgImage)

    avgRadIntegPlot = XYPlotData(goodEvents, "RADIAL INTEGRATION OF AVERAGE CSPAD", \
                                 radialPointsAvg,angularIntegAvg)
    socket.send("radintegavg", zmq.SNDMORE)
    socket.send_pyobj(avgRadIntegPlot)


