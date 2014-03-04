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

tt_hist = TTHistogram.load(sys.argv[1])

tt_hist.tau_hist()
tt_hist.plot(q_geom, 0.75, 3.0, 101, mask)


