
import sys
import numpy as np
from utilities import *

tt_hist = TTHistogram.load(sys.argv[1])
q_values = np.load('/reg/neh/home2/tjlane/analysis/xppb0114/averages/q_values.npy')
tt_hist.plot(q_values, 0.75, 3.5)

