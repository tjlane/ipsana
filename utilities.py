
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
            self.n_bins = n_bins
            self._bin_factor = float(self.n_bins-1) / self.q_values.max()
        else:
            self._bin_factor = 25.0
            self.n_bins = (self.q_values.max() * self._bin_factor) + 1
        
        self._bin_assignments = np.floor( q_values * self._bin_factor ).astype(np.int32)
        self._normalization_array = (np.bincount( self._bin_assignments.flatten(), weights=self.mask.flatten() ) \
                                    + 1e-100).astype(np.float)

        assert self.n_bins == self._bin_assignments.max() + 1
        self._normalization_array = self._normalization_array[:self.n_bins]
        
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
   
        assert bin_values.shape[0] == self.n_bins 
    
        return bin_values

    @property
    def bin_centers(self):
        return np.arange(self.n_bins) / self._bin_factor


class TTHistogram(object):

    def __init__(self, tau_min, tau_max, n_tau_bins):

        self._bin_cutoffs = np.linspace(tau_min, tau_max, n_tau_bins)
        self._n_shots = np.zeros(n_tau_bins-1, dtype=np.int)

        self._binned_intensities = np.zeros((n_tau_bins, 32, 185, 388))

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

         bin_index = ((self.bin_cutoffs >= tau-self.bin_width) * (self.bin_cutoffs < tau))[1:]
         if np.sum(bin_index) < 1:
             return
         elif np.sum(bin_index) > 1:
             print '%f in many bins?' % tau
             return

         assert data.shape == (32, 185, 388)
         print 'adding -->', tau
         self._n_shots[bin_index] += 1

         self._binned_intensities[bin_index] += data

         #self._binned_intensities[bin_index] *= (n-1)/float(n)
         #self._binned_intensities[bin_index] += (1.0/float(n)) * data

         return


    def histogram(self, q_geom, n_bins, mask=None):
        
        if mask == None:
            mask = np.ones((32, 185, 388), dtype=np.int)

        ra = RadialAverager(q_geom, mask, n_bins)
        
        histogram = np.zeros((len(self._n_shots), n_bins))

        for i in range(len(self._n_shots)):
            histogram[i,:] = ra(self._binned_intensities[i] / (float(self._n_shots[i]) + 1e-300))

        q_values = ra.bin_centers

        return q_values, histogram


    def plot(self, q_geom, q_min, q_max, n_q_bins, mask=None):

        q_values, hist_data = self.histogram(q_geom, n_q_bins, mask)

        inds = (q_values > q_min) * (q_values < q_max)
        qv = q_values[inds]

        x, y = np.meshgrid(qv, self.bin_centers)
        z = hist_data[:,inds]

        # compute relative change
        calib = np.copy(z[0,:][None,:])
        z -= calib
        #z /= calib

        fig = plt.figure()

        ax = fig.add_subplot(111, projection='3d')
        ax.plot_wireframe(x, y, z)

        ax.set_xlabel(r'Momentum Transfer $(\AA^{-1})$')
        ax.set_ylabel(r'Time Delay $\tau$ (ps)')
        ax.set_zlabel('Relative Intensity Change (%)')

        plt.show()

        return


    def timebin_dist_plot(self):
        plt.figure()
        plt.hist(self._n_shots)
        plt.show()
        return


    def tau_hist(self):
        plt.figure()
        plt.plot(self.bin_centers, self._n_shots)
        plt.xlabel(r'$\tau$ (ps)')
        plt.ylabel('Shots')
        plt.show()
        return


    def combine(self, other):
        if not isinstance(other, TTHistogram):
            raise TypeError('Can only combine with other TTHistogram objects')
        if not np.all(self._bin_cutoffs == other._bin_cutoffs):
            raise ValueError('bin cutoffs for both TT Histograms must match to'
                             ' combine them!')

        self._n_shots += other._n_shots
        self._binned_intensities += other._binned_intensities
        return


    def save(self, filename):
        io.saveh(filename,
                 bin_cutoffs = self._bin_cutoffs,
                 n_shots     = self._n_shots,
                 intensities = self._binned_intensities)
        print('Wrote %s to disk.' % filename)
        return


    @classmethod
    def load(cls, filename):
        hdf = io.loadh(filename)
        klass = cls(0.0, 1.0, 2)
        klass._bin_cutoffs = hdf['bin_cutoffs']
        klass._n_shots     = hdf['n_shots']
        klass._binned_intensities = hdf['intensities']
        return klass


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



