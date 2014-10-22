
import psana

class Event(object):
    """
    A container class that lightly wraps a psana event and exposes everything
    we want...
    """
    
    def __init__(self, psana_event):
        self.psana_event = psana_event
        self.ds1_src = psana.Source('DetInfo(XppGon.0:Cspad.0)')
        self.ds2_src = psana.Source('DetInfo(XppGon.0:Cspad.0)')
        self.evr_src = psana.Source('DetInfo(NoDetector.0:Evr.0)')
        return
        
        
    @property
    def pumpprobe_delay(self):
        return
        
    @property
    def pump_laser_status(self):
        return
        
    @property
    def xfel_status(self):
        return
        
    @property
    def pulse_energy_eV(self):
        return
        
    @property
    def ds1_distance_mm(self):
        return
        
    @property
    def ds2_distance_mm(self):
        return
        
    @property
    def ds1_instensities(self):
        return
    
    @property
    def ds2_instensities(self):
        return

