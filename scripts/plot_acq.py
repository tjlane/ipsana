
import time

from matplotlib import pyplot as plt
import psana

from pypad import cspad
from pypad.read import enforce_raw_img_shape as eris

plt.ion()

psana.setConfigFile('/reg/data/ana14/cxi/cxif7214/res/cfg/cxif7214.cfg')

cspad_ds1_src  = psana.Source('DetInfo(CxiDs1.0:Cspad.0)')
evr_src = psana.Source('DetInfo(NoDetector.0:Evr.0)')
aqr_src = psana.Source('DetInfo(CxiEndstation.0:Acqiris.0)')

ds = psana.DataSource('exp=cxif7214:run=64')

d = cspad.CSPad.load('/reg/data/ana14/cxi/cxif7214/scratch/averages/approx.cspad')


for evt in ds.events():

    a = evt.get(psana.Acqiris.DataDescV1, aqr_src)
    ch = a.data(5)

    fifos = evt.get(psana.EvrData.DataV3, evr_src).fifoEvents()
    evrs = [e.eventCode() for e in fifos]

    ds1 = evt.get(psana.ndarray_float32_3, cspad_ds1_src, 'calibrated_ndarr')


    im = plt.imshow( d(eris( ds1 )), vmin=0, vmax=2)
    if 183 in evrs:
        im.set_data( d(eris( ds1 )) )
        plt.draw()
        time.sleep(0.1)


