
import psana

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
                 ds1_offset=0.0, ds2_offset=0.0, energy_offset=0.0):
                 
        self.psana_event = psana_event
        
        self.ds1_src = psana.Source('DetInfo(CxiDs1.0:Cspad.0)')
        self.ds2_src = psana.Source('DetInfo(CxiDs2.0:Cspad.0)')
        self.evr_src = psana.Source('DetInfo(NoDetector.0:Evr.0)')
        
        self.ds1_offset = ds1_offset
        self.ds2_offset = ds2_offset
        self.energy_offset = energy_offset
        
        self._get_EVR_codes()

        return
        
    def _get_EVR_codes(self):
        fifos = self.psana_event.get(psana.EvrData.DataV3, self.evr_src).fifoEvents()
        self.evr_codes = [e.eventCode() for e in fifos]
        return
        
    @property
    def pumpprobe_delay_fs(self):
        return 0.0 
        
    @property
    def pump_laser_status(self):
        if (183 in self.evr_codes) and not (188 in self.evr_codes):
            status = 1
        elif not (183 in self.evr_codes) and (188 in self.evr_codes):
            status = 0
        else:
            #status = 'unknown'
            status = 0
        return status
        
    @property
    def xfel_status(self):

        if (187 in self.evr_codes) and not (189 in self.evr_codes):
            status = 1
        elif not (187 in self.evr_codes) and (189 in self.evr_codes):
            status = 0
        else:
            #status = 'unknown'
            status = 0 

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
    def ds1_intensities(self, corrected=True):
        if corrected:
            image = self.psana_event.get(psana.ndarray_float32_3, self.ds1_src, 
                                        'calibrated_ndarr')
        else:
            cspad = self.psana_event.get(psana.CsPad.DataV2, self.ds1_src)
            image = np.vstack([ cspad.quads(i).data() for i in range(4) ])
        if image == None: print 'warning: ds1 is NONE'
        return image
    
    @property
    def ds2_intensities(self, corrected=True):
        if corrected:
            image = self.psana_event.get(psana.ndarray_float32_3, self.ds2_src, 
                                         'calibrated_ndarr')
        else:
            cspad = self.psana_event.get(psana.CsPad.DataV2, self.ds2_src)
            image = np.vstack([ cspad.quads(i).data() for i in range(4) ])
        if image == None: print 'warning: ds2 is NONE'
        return image


