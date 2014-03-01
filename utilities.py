
import os
import cPickle

import numpy as np
from mdtraj import io

from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from pypad.read import enforce_raw_img_shape

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
                  162 : (0, 0)} # all dark
    
    
    # First experiment :: xppi0
    # uv on  & x-ray on: 140
    # uv off & x-ray on: 41, 67, 141
    # all dark: 162, 163
    # return is (xfel_status, uv_status)
    # evr_legend = {140 : (1, 1), # uv on/xray on
    #               41  : (1, 0), # uv off/xray on
    #               67  : (1, 0),
    #               141 : (1, 0),
    #               162 : (0, 0), # all dark
    #               163 : (0, 0)} 
    #               
    
    xfel_status, uv_status = (0,0) # default if no EVR code matches
    for fifoEvent in evr_data.fifoEvents():
        if fifoEvent.eventCode() in evr_legend.keys():
            xfel_status, uv_status = evr_legend[fifoEvent.eventCode()]
            #print "\tEVR EVENT: %d || XRAY: %d || UV: %d" % (fifoEvent.eventCode(), xfel_status, uv_status)

    return xfel_status, uv_status
    
    
class RadialAverager(object):
    
    def __init__(self, q_values, mask, n_bins=101):
        """
        Parameters
        ----------
        q_values : np.ndarray (float)
            For each pixel, this is the momentum transfer value of that pixel
        mask : np.ndarray (int)
            A boolean (int) saying if each pixel is masked or not
        n_bins : int
            The number of bins to employ. If `None` guesses a good value.
        """
        
        self.q_values = enforce_raw_img_shape(q_values)
        self.mask = enforce_raw_img_shape(mask).astype(np.int)
        self.n_bins = n_bins
        
        # figure out the number of bins to use
        if n_bins != None:
            self._bin_factor = float(self.n_bins-1) / self.q_values.max()
        else:
            self._bin_factor = 25.0
        
        self._bin_assignments = np.floor( q_values * self._bin_factor ).astype(np.int32)
        self._normalization_array = (np.bincount( self._bin_assignments.flatten(), weights=self.mask.flatten() ) \
                                    + 1e-100).astype(np.float)[:bin_values.shape[0]]
        
        return
    
    def __call__(self, image):
        """
        Bin pixel intensities by their momentum transfer.
        
        Parameters
        ----------            
        image : np.ndarray
            The intensity at each pixel, same shape as pixel_pos


        Returns
        -------
        bin_centers : ndarray, float
            The q center of each bin.

        bin_values : ndarray, int
            The average intensity in the bin.
        """

        image = enforce_raw_img_shape(image)
        
        if not (image.shape == self.q_values.shape):
            raise ValueError('`image` and `q_values` must have the same shape')
        if not (image.shape == self.mask.shape):
            raise ValueError('`image` and `mask` must have the same shape')

        weights = image.flatten() * self.mask.flatten()
        bin_values = np.bincount(self._bin_assignments.flatten(), weights=weights)
        bin_values /= self._normalization_array
    
        bin_centers = np.arange(bin_values.shape[0]) / self._bin_factor
    
        return bin_centers, bin_values


class TTHistogram(object):

    def __init__(self, tau_min, tau_max, n_tau_bins, n_q_bins):

        self._n_q_bins = n_q_bins
        self._bin_cutoffs = np.linspace(tau_min, tau_max, n_tau_bins)
        self._histogram = np.zeros(( n_tau_bins, n_q_bins ))
        self._n_shots = np.zeros(n_tau_bins, dtype=np.int)

        return

    @property
    def bin_cutoffs(self):
        return self._bin_cutoffs

    @property
    def bin_centers(self):
        return self.bin_cutoffs[1:] - self.bin_width/2.0

    @property
    def bin_width(self):
         return self.bin_cutoffs[1] - self.bin_cutoffs[0]

    def add_data(self, data, tau):

         if tau == 0.0:
             return # bad data

         bin_index = (self.bin_cutoffs >= tau-self.bin_width) * (self.bin_cutoffs < tau)
         if np.sum(bin_index) < 1:
             return
         elif np.sum(bin_index) > 1:
             print '%f in many bins?' % tau
             return

         assert len(data) == self._n_q_bins
         print 'adding -->', tau
         self._n_shots[bin_index] += 1
         n = int(self._n_shots[bin_index])
         self._histogram[bin_index] *= (n-1)/float(n)
         self._histogram[bin_index] += (1.0/float(n)) * data

         return


    def plot(self, q_values, q_min, q_max):

        inds = (q_values > q_min) * (q_values < q_max)
        qv = q_values[inds]

        x, y = np.meshgrid(qv, self.bin_cutoffs)
        z = self._histogram[:,inds]
        z -= z[0,:][None,:]

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        #ax.plot_wireframe(x, y, z)
        ax.plot_surface(x, y, z)
        plt.show()

        return


    def _to_serial(self):
        """ serialize the object to an array """
        s = np.array( cPickle.dumps(self) )
        s.shape=(1,) # a bit nasty...
        return s


    @classmethod
    def _from_serial(self, serialized):
        """ recover a Detector object from a serialized array """
        if serialized.shape == (1,):
            serialized = serialized[0]
        d = cPickle.loads( str(serialized) )
        return d


    def save(self, filename, overwrite=False):
        io.saveh(filename, detector=self._to_serial())
        print('Wrote %s to disk.' % filename)
        return


    @classmethod
    def load(cls, filename):
        hdf = io.loadh(filename)
        return cls._from_serial(hdf['detector'])


def normalize(q_values, intensities, q_min=0.5, q_max=3.5):
    assert q_values.shape == intensities.shape
    inds = (q_values > q_min) * (q_values < q_max)
    factor = float(np.sum(intensities[inds])) / float(np.sum(inds))
    return intensities / factor

def differential_integral(laser_on, laser_off, q_values, q_min=1.0, q_max=2.5):
    percent_diff = (laser_on - laser_off) / laser_on
    inds = (q_values > q_min) * (q_values < q_max)
    di = np.abs(np.sum( precent_diff[inds] ))
    return di



