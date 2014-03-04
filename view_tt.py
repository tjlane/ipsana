#!/usr/bin/env python

"""
Get a bunch of TT Histograms, combine them, and make beautiful plots :)
"""

import os
import sys
from glob import glob

import numpy as np

from utilities import *


q_geom = np.load('/reg/neh/home2/tjlane/analysis/xppb0114/geometries/v2/q_geom.npy')
mask = np.load('/reg/neh/home2/tjlane/analysis/xppb0114/geometries/v2/mask_v2.npy')

print 'Searching: %s' % sys.argv[1]
hist_filenames = glob(os.path.join(sys.argv[1], '*.h5'))
n_files = len(hist_filenames)
print 'Found: %d files...' % n_files

tt_hist = TTHistogram.load(hist_filenames[0])

for i,filename in enumerate(hist_filenames[1:]):
    print '%d/%d -- %s' % (i+1, n_files, filename)
    tt_hist.combine( TTHistogram.load(filename) )


tt_hist.tau_hist()
tt_hist.plot(q_geom, 0.75, 3.0, 101, mask)


