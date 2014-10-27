
"""
"""

import numpy as np
import psana


global ds1_src
ds1_src = psana.Source('DetInfo(CxiDs1.0:Cspad.0)')
global ds2_src
ds2_src = psana.Source('DetInfo(CxiDs2.0:Cspad.0)')
global evr_src
evr_src = psana.Source('DetInfo(NoDetector.0:Evr.0)')
global dark
#dark = np.load('/reg/d/psdm/cxi/cxif7214/scratch/averages/run107_dark.npy')
dark = np.load('/reg/neh/home2/tjlane/opt/minitti_beamtime/run107_dark.npy')
dark += np.load('/reg/neh/home2/tjlane/opt/minitti_beamtime/run101_ds1_average.npy')


def print_evrs(run_num):
    ds = psana.DataSource('exp=cxif7214:run=39')
    for evt in ds.events(): 
        print [e.eventCode() for e in evt.get(psana.EvrData.DataV3, psana.Source('DetInfo(NoDetector.0:Evr.0)')).fifoEvents()]
    return


class Event(object):
    """
    A container class that lightly wraps a psana event and exposes everything
    we want...
    
    EVR Codes
    ---------
    183 : Laser On
    184 : Mol Beam On
    185 : CXI Pulse Picker
    187 : X-ray beam on
    188 : Laser Off
    189 : X-ray beam off
    
    """
    
    def __init__(self, psana_event, epics, 
                 ds1_offset=0.0, ds2_offset=0.0, energy_offset=0.0,
                 corrected=True):
                 
        self.psana_event = psana_event
        self.corrected = corrected        

        self.ds1_offset = ds1_offset
        self.ds2_offset = ds2_offset
        self.energy_offset = energy_offset
        
        self._get_EVR_codes()

        return
        
    def _get_EVR_codes(self):
        if self.psana_event:
            fifos = self.psana_event.get(psana.EvrData.DataV3, evr_src).fifoEvents()
            self.evr_codes = [e.eventCode() for e in fifos]
        else:
            self.evr_codes = []
        return
        
    @property
    def pumpprobe_delay_fs(self):
        return 0.0 
        
    @property
    def pump_laser_status(self):
        if (183 in self.evr_codes) and not (184 in self.evr_codes):
            status = 1
        elif not (183 in self.evr_codes) and (184 in self.evr_codes):
            status = 0
        else:
            status = 'unknown'
        return status
        
    @property
    def xfel_status(self):
        if (187 in self.evr_codes) and not (189 in self.evr_codes):
            status = 1
        elif not (187 in self.evr_codes) and (189 in self.evr_codes):
            status = 0
        else:
            status = 'unknown'
        return status
        
    @property
    def pulse_energy_eV(self):
        E = self.epics.value('SIOC:SYS0:ML00:AO627')
        return E + self.energy_offset
        
    @property
    def ds1_distance_mm(self):
        z = self.epics.value('CXI:DS1:MMS:06.RBV')
        return z + self.ds2_offset
        
    @property
    def ds2_distance_mm(self):
        z = self.epics.value('CXI:DS2:MMS:07.RBV')
        return z + self.ds2_offset
        
    @property
    def polarization_status(self):
        return 0
        
    @property
    def ds1_intensities(self):
        if self.corrected:
            image = self.psana_event.get(psana.ndarray_float32_3, ds1_src, 
                                        'calibrated_ndarr')
        else:
            cspad = self.psana_event.get(psana.CsPad.DataV2, ds1_src)
            if cspad:
                image = np.vstack([ cspad.quads(i).data() for i in range(4) ])
                image -= dark
            else:
                image = None
        if image == None:
            print 'warning: ds1 is NONE'
        return image
    
    @property
    def ds2_intensities(self):
        if self.corrected:
            image = self.psana_event.get(psana.ndarray_float32_3, ds2_src, 
                                         'calibrated_ndarr')
        else:
            cspad = self.psana_event.get(psana.CsPad.DataV2, ds2_src)
            if cspad:
                image = np.vstack([ cspad.quads(i).data() for i in range(4) ])
            else:
                image = None
        if image == None: print 'warning: ds2 is NONE'
        return image


