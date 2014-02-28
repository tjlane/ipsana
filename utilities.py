
import numpy as np

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
    
    def __init__(self, q_values, mask, n_bins=100):
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
        
        self.q_values = q_values
        self.mask = mask
        self.n_bins = n_bins
        
        # figure out the number of bins to use
        if n_bins != None:
            bin_factor = float(self.n_bins) / self.q_values.max()
        else:
            bin_factor = 25.0
        
        self._bin_assignments = np.floor( q_values * bin_factor ).astype(np.int32)
        
        self.mask = enforce_raw_img_shape(self.mask) 
        self.mask = self.mask.astype(np.int).flatten()
        
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
        
        if not (image.shape == self.q_values.shape):
            raise ValueError('`image` and `q_values` must have the same shape')
        if not (image.shape == self.mask.shape):
            raise ValueError('`image` and `mask` must have the same shape')
    
        image = enforce_raw_img_shape(image)

        weights = image.flatten() * self.mask
        bin_values = np.bincount(self._bin_assignments.flatten(), weights=weights)
        bin_values /= (np.bincount( self._bin_assignments.flatten(), weights=self.mask ) \
                                    + 1e-100).astype(np.float)[:bin_values.shape[0]]
    
        bin_centers = np.arange(bin_values.shape[0]) / bin_factor
    
        return bin_centers, bin_values


def normalize(q_values, intensities, q_min=0.5, q_max=3.5):
    assert q_values.shape == intensities.shape
    inds = (q_values > q_min) * (q_values < q_max)
    factor = float(np.sum(intensities[inds])) / float(np.sum(inds))
    return intensities / factor





