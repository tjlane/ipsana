
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
        
        return
        
    def _get_EVR_codes(self):
        # evt.get(psana.EvrData.DataV3, evr_src)
        # self.evr_codes = 
        return
        
    @property
    def pumpprobe_delay_fs(self):
        return
        
    @property
    def pump_laser_status(self):
        return
        
    @property
    def xfel_status(self):
        return
        
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
        return
        
    @property
    def ds1_instensities(self, corrected=True):
        if corrected:
            image = self.psana_event.get(psana.ndarray_float32_3, self.ds1_src, 
                                        'calibrated_ndarr')
        else:
            cspad = self.psana_event.get(psana.CsPad.DataV2, self.ds1_src)
            image = np.vstack([ cspad.quads(i).data() for i in range(4) ])
        return image
    
    @property
    def ds2_instensities(self, corrected=True):
        if corrected:
            image = self.psana_event.get(psana.ndarray_float32_3, self.ds2_src, 
                                         'calibrated_ndarr')
        else:
            cspad = self.psana_event.get(psana.CsPad.DataV2, self.ds2_src)
            image = np.vstack([ cspad.quads(i).data() for i in range(4) ])
        return image

