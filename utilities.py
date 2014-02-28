
import numpy as np

# Use first EVR code
# uv on  & x-ray on: 140
# uv off & x-ray on: 41, 67, 141
# all dark: 162, 163


def update_average(n, A, B):
    '''
    updates a numpy matrix A that represents an average over the previous n-1 shots
    by including B into the average, B being the nth shot
    '''
    A *= (n-1)/float(n)
    A += (1.0/float(n))*B
    return

    
def pumpprobe_status(evr_data):
    """
    """
    
    # Second experiment :: 
    # uv on  & x-ray on : 90
    # no UV : 91
    # no xfel : 62
    # return is (xfel_status, uv_status)
    evr_legend = {90  : (1, 1), # uv on/xray on
                  91  : (1, 0), # uv off/xray on
                  162 : (0, 0)] # all dark
    
    
    # First experiment :: xppi0
    # uv on  & x-ray on: 140
    # uv off & x-ray on: 41, 67, 141
    # all dark: 162, 163
    # return is (xfel_status, uv_status)
    evr_legend = {140 : (1, 1), # uv on/xray on
                  41  : (1, 0), # uv off/xray on
                  67  : (1, 0),
                  141 : (1, 0),
                  162 : (0, 0), # all dark
                  163 : (0, 0)} 
                  

    xfel_status, uv_status = (None, None)
    for fifoEvent in evr_data.fifoEvents():
        if fifoEvent.eventCode() in evr_legend.keys():
            xfel_status, uv_status = evr_legend[fifoEvent.eventCode()]
            print "\tEVR EVENT: %d || XRAY: %d || UV: %d" % (fifoEvent.eventCode(), xfel_status, uv_status)


    return xfel_status, uv_status