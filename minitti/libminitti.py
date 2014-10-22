#!/usr/bin/env python

"""
A grab-bag of algorithms we'll need...

NOTES
-----

--> cspad data
Standard data shape is (2, 32, 185, 388), the first "2" is for DS1/DS2.

"""

import numpy as np


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
        bin_values /= self._normalization_array

        assert bin_values.shape[0] == self.n_bins

        return bin_values
    

    @property
    def bin_centers(self):
        return np.arange(self.n_bins) / self._bin_factor
        
        

class _TTHistogramBin(object):
    """
    This is a container class used by TTHistogram
    """
    
    
    def __init__(bin, num_laseroff_shots=0, num_laseron_shots=0,
                 laseron_avg=np.zeros(2, 32, 185, 388),
                 laseroff_avg=np.zeros(2, 32, 185, 388),
                 cspad_SminusP=np.zeros(2, 32, 185, 388)):
                 
        # set all values internally
        self.bin = bin
        self.num_laseroff_shots = num_laseroff_shots
        self.num_laseron_shots = num_laseron_shots
        
        self.laseron_avg = laseron_avg
        self.laseroff_avg = laseroff_avg
        
        self.cspad_SminusP = cspad_SminusP
        
        return
        
        
    @property
    def num_shots(self):
        return self.num_laseroff_shots + self.num_laseron_shots
    
        
    def add_shot(cspad_intensities, laser_status, polarization):
        
        # laser status
        if laser_status == 0:
            self.num_laseroff_shots += 1
            update_average(self.num_laseroff_shots, self.laseroff_avg, 
                           cspad_intensities)
                           
        elif laser_status == 1:
            self.num_laseron_shots += 1
            update_average(self.num_laseron_shots, self.laseron_avg, 
                           cspad_intensities)
                           
        else:
            raise ValueError('laser_status invalid value: %d' % laser_status)
            
            
        # polarization
        if polarization not in [-1, 1]:
            raise ValueError('polarization not in [-1, 1]')
        self.cspad_SminusP += polarization * cspad_intensities
        
        return
    

class TTHistogram(object):
    
    def __init__(self, ds1_geometry, ds2_geometry, bin_width=5.0):
        """
        
        Parameters
        ----------
        geometry : np.ndarray
            2D Numpy array, (32*185*388) x 3, describing the reciprocal space
            coordinates of each pixel in the CSPAD detector (q, theta, phi).
            See thor.Detector.recpolar.
        
        bin_width : float
            The width of each bin, in femtoseconds. Default is 5fs.
        """
        
        self._bin_width = bin_width
        
        for g in [ds1_geometry, ds2_geometry]:
            if not g.shape == (32*185*388, 3):
                raise ValueError('`geometry` is not the correct shape!')
        
        self._ds1_geometry = ds1_geometry
        self._ds2_geometry = ds2_geometry
        
        # we're going to keep track of the bins using an "index_multiplier",
        # which is a unique index for each bin. The index_multiplier * bin_width
        # gives the left edge of the bin (inclusive). For example, if the
        # bin_width = 5.0, and the value tau=47, then that tau goes in the bin
        # spanning 45.0-50.0, who's index_multiplier is 45/5 = 9
        
        self._tt_hist_bins = {}

        return
        
        
    @property
    def ds1_geometry(self):
        return self._ds1_geometry
    

    @property
    def ds2_geometry(self):
        return self._ds2_geometry
    
        
    @property
    def bin_width(self):
        return self._bin_width
    
    @property
    def _index_multipliers_active(self):
        return np.array( sorted(self._tt_hist_bins.keys()) )
    

    @property
    def bin_left_edges(self):
        return self._index_multipliers_active.astype(np.float) * self.bin_width
    

    @property
    def bin_centers(self):
        return self.bin_left_edges + self.bin_width / 2.0
    

    @property
    def num_shots_per_bin(self):
        shots_per_bin = []
        for k in self._index_multipliers_active()
            shots_per_bin.append(self._tt_hist_bins[k])
        return np.array(shots_per_bin)
    
        
    def _determine_index_multipier(tau):
        return int( floor(tau / self._bin_width) )
    

    def _initialize_new_hist_bin(self, index_multiplier, **kwargs):
        
        if index_multiplier in self._index_multipliers_active:
            raise RuntimeError('trying to add existing bin')
            
        self._tt_hist_bins[index_multiplier] = TTHistogramBin(index_multiplier, 
                                                              **kwargs)
        
        return
    

    def add_image(self, tau, laser_status, polarization, cspad_intensities):
        """
        Add a tagged CSPAD image
        """
        
        im = self._determine_index_multipier(tau)
        
        if not im in self._index_multipliers_active:
            self._initialize_new_hist_bin(im)
            
        # now update the appropriate bin
        self._tt_hist_bins[im].add_shot(cspad_intensities, 
                                        laser_status, 
                                        polarization)
        
            
        
        return
        
        
    def coarse_grain(new_bin_width):
        
        if new_bin_width <= self.bin_width:
            raise ValueError('new bin width is smaller or equal to old')
            
        raise NotImplementedError()
        
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
        #io.saveh(filename,
        #         bin_cutoffs = self._bin_cutoffs,
        #         n_shots     = self._n_shots,
        #         intensities = self._binned_intensities)
        
        f = h5py.File(filename)
        f['/bin_cutoffs'] = self._bin_cutoffs
        f['/n_shots']     = self._n_shots
        f['/intensities'] = self._binned_intensities
        f.close()        

        return


    @classmethod
    def load(cls, filename):
        return klass
        
        
def update_average(n, A, B):
    """
    updates a numpy matrix A that represents an average over the previous n-1 shots
    by including B into the average, B being the nth shot
    """
    if n == 0:
        A += B
        return
    else:
        A *= (n-1)/float(n)
        A += (1.0/float(n))*B
        return
        
        
def normalize(q_values, intensities, q_min=0.5, q_max=3.5):
    """
    Crop and normalize an I(q) vector st. the area under the curve is one.
    """
    assert q_values.shape == intensities.shape
    inds = (q_values > q_min) * (q_values < q_max)
    factor = float(np.sum(intensities[inds])) / float(np.sum(inds))
    return intensities / factor


def differential_integral(laser_on, laser_off, q_values, q_min=1.0, q_max=2.5):
    """
    Compute the following
    
        DI = \int_{q0}^{q1} | lazer_on(q) - lazer_off(q) | dq
        
    Useful for tracking changes in the scattering.
    """
    percent_diff = (laser_on - laser_off) / laser_on
    inds = (q_values > q_min) * (q_values < q_max)
    di = np.abs(np.sum( precent_diff[inds] ))
    return di
    
