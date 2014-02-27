#!/usr/bin/env python

# ------------------------------------
# Setup PYTHON path to add latest ZMQ
# and Dan D's path
import os
import sys
import time
# ------------------------------------

import numpy as np
import psana 

# ZMQ RELATED MODULES
import zmq
from psdata import ImageData
from psdata import XYPlotData

####################
# ZMQ SETUP 
context = zmq.Context()
socket = context.socket(zmq.PUB)
socket.setsockopt(zmq.SNDHWM, 10)
socket.bind("tcp://*:%d" % 12322)
####################

# ------- CONFIG -------
# PSANA CONFIG FILE - this must be set before a datasource is created
print os.path.abspath(__file__)
basedir = os.path.split(os.path.abspath( __file__ ))[0]
configFileName = os.path.join(basedir,"minitti.cfg")
assert os.path.exists(configFileName), "config file not found, looked for: %s" % configFileName
psana.setConfigFile(configFileName)
psana.setOption('psana.l3t-accept-only',0)

imgkey = 'image0'  

# Define experiment, run. Shared memory syntax (for SXR): shmem=0_1_SXR.0:stop=no
ds = psana.DataSource('exp=xppi0613:run=275') 

cspadSrc  = psana.Source('DetInfo(XppGon.0:Cspad.0)')
evrSrc    = psana.Source('DetInfo(NoDetector.0:Evr.0)')
FEEGasSrc = psana.Source('BldInfo(FEEGasDetEnergy)')
statusPrintout = 100      # print a status every on the Nth shot
updatePlots = 100         # carry out radial integration and update plots every N good shots
laserOnEvrCode = 41      # for positive N, means laser is present when evr code is present
                         # set to -N if laser is on when evr code N is not present
noBeamEvrCode = 162
RadialBinPixelWidth = 5  # about how many pixels you want in a bin for radial integration
                         # 100 bins is 12.4 pixels.

rStart = 100.0             # where to start radial integration
rEnd =   400.0             # max for cspad is 1224 

ttpvs = ['TTSPEC:AMPL', 'TTSPEC:AMPLNXT', 'TTSPEC:FLTPOS', 
         'TTSPEC:FLTPOSFWHM', 'TTSPEC:FLTPOS_PS', 'TTSPEC:REFAMPL']

checkedThatPvsArePresent = False
ttvals = {}
ttbounds={}
ttbounds['laser_on'] = {}
ttbounds['laser_off'] = {}

'''
ttbins: {'laser_off': {'TTSPEC:REFAMPL': {'max': 4655.08984375,     'min': 2837.889892578125}, 
          'laser_on': {'TTSPEC:REFAMPL': {'max': 4655.08984375, 'min': 2837.889892578125}, 
                    'TTSPEC:FLTPOSFWHM': {'max': 2453.291267802251, 'min': 31.097311990086677}, 
                       'TTSPEC:AMPLNXT': {'max': 0.28004136113793265, 'min': 0.10136527900328929}, 
                     'TTSPEC:FLTPOS_PS': {'max': 0.7994734047079562, 'min': -0.6075994956608504}, 
                        'TTSPEC:FLTPOS': {'max': 944.9017135592403, 'min': 255.16009573139402}, 
                          'TTSPEC:AMPL': {'max': 1.9576389043495226, 'min': 0.1453040483727932}}, 
                    'TTSPEC:FLTPOSFWHM': {'max': 2453.291267802251, 'min': 31.097311990086677}, 
                       'TTSPEC:AMPLNXT': {'max': 0.28004136113793265, 'min': 0.10136527900328929}, 
                     'TTSPEC:FLTPOS_PS': {'max': 0.7994734047079562, 'min': -0.6075994956608504}, 
                        'TTSPEC:FLTPOS': {'max': 944.9017135592403, 'min': 255.16009573139402}, 
                          'TTSPEC:AMPL': {'max': 1.9576389043495226, 'min': 0.1453040483727932}}}
'''

for ky in ttbounds:
    for pv in ttpvs:
        ttbounds[ky][pv]={'min':9e32,'max':-9e32}
        ttvals[pv]=None

# ------ END CONFIG ---------

# ----- Helper functions
def updateAverage(n,A,B):
    '''updates a numpy matrix A that represents an average over the previous n-1 shots
    by including B into the average, B being the nth shot
    '''
    A *= (n-1)/float(n)
    A += (1.0/float(n))*B

def getBeamLaserFromEvr(evrData,laserOnEvrCode,noBeamEvrCode):
    '''returns beamPresent and laserOn
    ARGS:
      evrData        - as returned from psana event get function
      laserOnEvrCode - if positive, evr code that means laser is present, if negative,
                       means laser is not present when evr code is present
      noBeamEvrCode  - evr code present when no beam fires
    '''
    evrData = evt.get(psana.EvrData.DataV3, evrSrc)
    if evrData is None:  # evr will always be present, this shouldn't be 
        return None,None # neccessary
    
    beamPresent = True
    laserCodePresent=False
    for fifoEvent in evrData.fifoEvents():
        if fifoEvent.eventCode() == abs(laserOnEvrCode):
            laserCodePresent = True
        elif fifoEvent.eventCode() == noBeamEvrCode:
            beamPresent = False

    if (laserCodePresent and laserOnEvrCode>0) or (not laserCodePresent and laserOnEvrCode<0):
        laserOn = True
    else:
        laserOn = False

    return beamPresent, laserOn # this is XFEL/UV

# ------ end helper functions

eventCounter = 0
goodEvents = 0 
laserEvents = 0
noLaserEvents = 0
badEvents = 0   
avgCspadAll = None
avgCspadLaser = None
avgCspadNoLaser = None

radIntegrator = AngularIntegrator()
epics = ds.env().epicsStore()
calib = ds.env().calibStore()

onoffDict = {True:'on',False:'off'}
for evt in ds.events():
    t0 = time.time()
    eventCounter += 1

    cspadCalib = evt.get(psana.ndarray_float32_2, cspadSrc, imgkey)
    
    if not checkedThatPvsArePresent:
        for pvName in ttpvs:
            assert pvName in epics.pvNames(), "time tool pv %s not found in epics store" % pvName
        checkedThatPvsArePresent = True

    for pvName in ttpvs:
        ttvals[pvName] = epics.value(pvName)        

    evrData = evt.get(psana.EvrData.DataV3, evrSrc)
    if evrData is None:
        badEvents += 1
        continue

    beamPresent, laserOn = getBeamLaserFromEvr(evrData,laserOnEvrCode,noBeamEvrCode)

    onoff = onoffDict[laserOn]
    for pvName in ttpvs:
        ttbounds['laser_%s'%onoff][pvName]['min'] = min(ttbounds['laser_on'][pvName]['min'], ttvals[pvName])
        ttbounds['laser_%s'%onoff][pvName]['max'] = max(ttbounds['laser_on'][pvName]['max'], ttvals[pvName])

    FEEGasValues = None
    FEEGasMsg = " no FEEGas this event"
    FEEGas = evt.get(psana.Bld.BldDataFEEGasDetEnergy, FEEGasSrc)
    if FEEGas is not None:
        FEEGasValues = {'11':FEEGas.f_11_ENRC(), 
                        '21':FEEGas.f_12_ENRC(), 
                        '12':FEEGas.f_21_ENRC(), 
                        '22':FEEGas.f_22_ENRC() }
        FEEGasMsg = "FEEGas: 11: %f" % FEEGasValues['11']

    if (eventCounter % statusPrintout == 0) :
        print "cspad_radial - Total good events: ", goodEvents, " bad: ", badEvents
        print "cspad_radial - %s" % FEEGasMsg

    if cspadCalib is None:
        badEvents +=1
        continue

    goodEvents += 1

    if avgCspadAll is None:
        avgCspadAll = np.zeros(cspadCalib.shape)
        avgCspadLaser = np.zeros(cspadCalib.shape)
        avgCspadNoLaser = np.zeros(cspadCalib.shape)
        cspadCenterX = int(cspadCalib.shape[0]/2)
        cspadCenterY = int(cspadCalib.shape[1]/2)
        radialBins = int(np.max(np.array(cspadCalib.shape[0])/2.0)/RadialBinPixelWidth)
        cspadMask = calib.get(psana.ndarray_int16_2, cspadSrc)
        #assert cspadMask is not None, "Did not get cspadMask from store. You need to run from release dir with CsPadPixCoords tag > V00-03-17"
        
#        import IPython
#        IPython.embed()

    updateAverage(goodEvents, avgCspadAll, cspadCalib)

    if laserOn and beamPresent:
        laserEvents += 1
        updateAverage(laserEvents, avgCspadLaser, cspadCalib)

    if not laserOn and beamPresent:
        noLaserEvents += 1
        updateAverage(noLaserEvents, avgCspadNoLaser, cspadCalib)
        
    if goodEvents % updatePlots != 0:
        continue

    # Perform angular integration
    radialPointsAvg,angularIntegAvg = radIntegrator.angularIntegration(avgCspadAll, 
                                      cspadCenterX, cspadCenterX, nbins=radialBins, 
                                      rStart=rStart, rEnd=rEnd, mask=cspadMask)

    radialPointsEvent,angularIntegEvent = radIntegrator.angularIntegration(cspadCalib, 
                                          cspadCenterX, cspadCenterX, nbins=radialBins, 
                                          rStart=rStart, rEnd=rEnd,  mask=cspadMask)

    radialPointsLaser,angularIntegLaser = radIntegrator.angularIntegration(avgCspadLaser, 
                                          cspadCenterX, cspadCenterX, nbins=radialBins, 
                                          rStart=rStart, rEnd=rEnd, mask=cspadMask)

    radialPointsNoLaser,angularIntegNoLaser = radIntegrator.angularIntegration(avgCspadNoLaser, 
                                              cspadCenterX, cspadCenterX, nbins=radialBins,  
                                              rStart=rStart, rEnd=rEnd, mask=cspadMask)

    angularIntegLaserNorm = angularIntegLaser / np.sum(angularIntegLaser)
    angularIntegNoLaserNorm = angularIntegNoLaser / np.sum(angularIntegNoLaser)
    uvDiff = (angularIntegLaserNorm - angularIntegNoLaserNorm)/(1.0 + np.abs(angularIntegLaserNorm) + np.abs(angularIntegNoLaserNorm))

    print "cspad_radial ttbounds " , ttbounds
    print "cspad_radial ttvals  ", ttvals

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

    evtRadIntegPlot = XYPlotData(goodEvents, "RADIAL INTEGRATION OF CURRENT CSPAD", \
                                 radialPointsEvent,angularIntegEvent)
    socket.send("radintegevt", zmq.SNDMORE)
    socket.send_pyobj(evtRadIntegPlot)

    laserRadIntegPlot = XYPlotData(laserEvents, "RADIAL INTEGRATION OF AVERAGE LASER CSPAD", \
                                 radialPointsLaser,angularIntegLaser)
    socket.send("radinteglaser", zmq.SNDMORE)
    socket.send_pyobj(laserRadIntegPlot)


    noLaserRadIntegPlot = XYPlotData(noLaserEvents, "RADIAL INTEGRATION OF AVERAGE NO LASER CSPAD", \
                                 radialPointsNoLaser,angularIntegNoLaser)
    socket.send("radintegnolaser", zmq.SNDMORE)
    socket.send_pyobj(noLaserRadIntegPlot)

    uvDiffPlot = XYPlotData(goodEvents, "(UVON-UVOFF)/(1+|UVON|+|UVOFF|)", \
                                 radialPointsLaser,uvDiff)
    socket.send("uvdiff", zmq.SNDMORE)
    socket.send_pyobj(uvDiffPlot)

    print "cspad_radial - time for analyzing event: %.4e" % (time.time()-t0,)

