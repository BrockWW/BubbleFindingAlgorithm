'''File containing the main class for the bubble finding algorithm. Each object created with the class
completely runs the bubble finding algorithm with no need to call additional methods. Any methods run
with the initalization of the object should not be ran again. Additional methods are included to calculate
different properties of the found bubbles, but all data on the image is stored in self variables so
additional claculations can be done outside of the class methods.'''

# imports
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
import numpy as np
from noise_algorithm import *
from scipy import ndimage
from skimage.measure import label, regionprops

class Bubble_Alg():
    def __init__(self, im, im_size, percentile, beam_size = 1, pixel_size = 1):
        '''Class to combine all bubble algorithm functions and output into one coherent object
        to increase ease of use and speed of analysis.
        
        Variables:
            - im: (2D array) galaxy column density image to be analyzed, needs to be square (H*cm^-2)
            - im_size: (float) physical size of one side of input image (kpc)
            - percentile: (integer) volume percentile used to find the bubbles in galactic disk
            - beam_size: (float) major axis of beam used to smooth image (or used in retrieving 
                         observation image) (arcsec)
            - pixel_size: (float) angular size of each pixel (arcsec/pixel_length)
            '''
        
        # checking if image is square
        if(im.shape[0] != im.shape[1]):
            return "Image must be square with same pixel count on each side."
        
        # defining self variables fo use in following methods
        self.im = im
        self.im_size = im_size
        self.pixel_area = (im_size/im.shape[0])**2
        self.perc = percentile
        self.beam_size = beam_size
        self.pixel_size = pixel_size
        
        # running methods to find bubbles
        self.remove_exterior()
        self.ext_bubble_algorithm()
        self.ext_find_bubbles()
        
    def remove_exterior(self):
        '''Function that smooths imput image to find exterior region for removal.

        Variables:
            - self: calling self variables for function

        Returns:
            - self.ext_im: input image with all exterior values set to np.nan
            '''

        # selecting size parameter for uniform_filter based on image data
        l = self.im.shape[0]
        weight_im = self.im*((self.im_size*1000*3.086e18/l)**2)*(1.67e-27)   # H/cm^2 -> kg
        stat = np.average(self.im, weights = weight_im)
        v = np.sum(self.im > stat)
        # finding appoximate radius of galaxy in pixels
        s = np.sqrt(v)/2

        # using scipy to smooth image into exterior and galaxy regions only
        filt_im = ndimage.uniform_filter(self.im, size = s)
        thresh_im = filt_im < np.average(filt_im)

        # labeling image components and returning image of exterior removed galaxy
        labeled_im = label(thresh_im)
        ext_im = np.where(labeled_im == 1, np.nan, 0)
        ext_rem_im = np.where(np.isnan(ext_im), np.nan, self.im)

        # assigning exterior removed image to class properties
        self.ext_rem_im = ext_rem_im

    def ext_bubble_algorithm(self):
        '''New bubble finding algorithm that utilizes skimage rather than scipy to try
        and improve the ability of the algorithm to deal with noise and correctly find
        bubbles.

        Variables:
            - self: calling self variables for function

        Returns:
            - self.labeled_im: image with all seperated regions labeled uniquely
            - self.info_dict: dictionary that contains significant amout of data about each 
              labeled region (area, label value, centroid, etc.)
            - self.thresh_used: column density threshold (H*cm^-2)
        '''

        # finding thresholds for bubbles
        colden_thresh = np.nanpercentile(self.ext_rem_im, self.perc)
        thresh_im = self.ext_rem_im < colden_thresh

        # seperating each bubble by labeling
        labeled_im = label(thresh_im)

        # reading information about labeled regions
        info_dict = regionprops(labeled_im)
        
        # assigning new parameters to the class
        self.labeled_im = labeled_im
        self.info_dict = info_dict
        self.colden_thresh = colden_thresh

    def ext_find_bubbles(self):
        '''Bubble classifying algoithm that uses the volume percentile threshold method.

        Variables:
            - self: calling self variables for function

        Returns:
            - bubble_arr: array containing the region label values classified as bubbles
            - nb_arr: array containing the region label values not chosen as bubbles
        '''

        # finding minimum limit for bubble sorting based on beam size
        min_size = np.sqrt(self.pixel_area) * 2*self.beam_size / self.pixel_size   # kpc

        self.lr = []
        self.labels = []
        for i in range(len(self.info_dict)):
            self.labels.append(self.info_dict[i]["label"])
            diam = 2*np.sqrt(self.info_dict[i]["area"]*self.pixel_area/np.pi)
            if(diam >= min_size):
                self.lr.append(self.info_dict[i]["label"])
                
        binary_im = np.where(np.isnan(self.ext_rem_im), 1, 0)
        border_im = np.ones_like(binary_im)

        # get the size of the input array in both dimensions
        x_size = np.shape(binary_im)[0]
        y_size = np.shape(binary_im)[1]

        # iterate through each value
        for i in range(x_size):
            for j in range(y_size):
                # select exterior
                if binary_im[i,j] == 0:
                    # get slice of all adjacent pixels
                    slice = binary_im[i-1:i+1,j-1:j+1]
                    # check if any in galaxy
                    if np.any(slice == 1):
                        # update borders array
                        border_im[i,j] = 0
        
        bubbles = []
        # looping through each found region, determining if on border and removing
        for val in self.lr:
            binary_bub_im = np.where(self.labeled_im == val, 1, 0)
            compare_bub_im = binary_bub_im*border_im
               
            # testing if bubble is on border or not
            if(np.all(compare_bub_im == binary_bub_im)):
                bubbles.append(val)

        # finding non-bubble labels
        not_bubbles = [nb for nb in self.labels if nb not in bubbles]
        
        # saving bubble and non-bubble region labels to class parameters
        self.bubble_arr = np.array(bubbles)
        self.nb_arr = np.array(not_bubbles)

    def ext_find_bubble_radius(self):
        '''Calculating found bubble approximate radius by assuming all bubbles are 
        perfectly circular.

        Variables:
            - self: calling self variables for function

        Returns:
            - self.radius_arr: array of calculated found bubble approximate radii (kpc)
        '''

        # initializing bubble diameter array
        radius_arr = np.zeros_like(self.bubble_arr)

        # looping through all bubbles selected from ext_find_bubbles
        for i in range(len(self.bubble_arr)):
            # correcting for label-index difference
            index = self.bubble_arr[i] - 1
            # calculating bubble approximate radius
            radius_arr[i] = np.sqrt(self.info_dict[index]["area"]*self.pixel_area/np.pi)

        # saving radius of bubble to class parameters
        self.radius_arr = radius_arr
        
    def ext_bubble_galactic_radius(self, center):
        '''Algorithm for calculating bubble galactic radius.

        Variables:
            - self: calling self variables for function
            - center: array denoting the central pixel location of the image [px, py]

        Returns:
            - self.bubble_gal_rad: array of calculated found bubble center distances from 
                                   galctic center (kpc)
        '''
        
        # ensuring input center coordinate is a numpy array for later calculation
        if(not isinstance(center, np.ndarray)):
            center = np.array(center)

        # initializing bubble galactic radius array
        self.bubble_gal_rad = np.zeros_like(self.bubble_arr, dtype = np.float64)
        
        # looping through all found bubbles and calculating galactic radius
        for i in range(len(self.bubble_arr)):
            index = self.bubble_arr[i] - 1
            bub_center = self.info_dict[index]["centroid"]
            
            pixel_rad = np.sqrt(np.sum((center - bub_center)**2))
            kpc_rad = pixel_rad*np.sqrt(self.pixel_area)
            self.bubble_gal_rad[i] = kpc_rad