'''Script that includes the function to plot the found bubbles after using the bubble finding
algorithm. Output plot will have no labels or title, and will require these to be included 
outside of the function by the user.'''

# imports 
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, ListedColormap
import numpy as np

def ext_bubble_plotting(alg_obj):
        '''Algorithm for plotting found bubbles from algorithm.

        Variables:
            - alg_object: object returned from calling bubble algorithm class

        Returns:
            - ax: image of bubbles found by algorithm with denoted exterior and 
                  galactic disk regions free of labels
            - bub_gal_fraction: fractional coverage of bubbles on galactic disk
        '''
        
        # masking out realistic bubbles sizes only to plot
        bub_im = alg_obj.labeled_im
        
        # selecting if we want to show all regions or just bubbles
        for num in alg_obj.bubble_arr:
            bub_im = np.where(bub_im == num, 1000, bub_im)
        bub_im = np.where(bub_im == 1000, 10, 0)
          
        # setting exterior region to unique value for plotting
        bub_im = np.where(np.isnan(alg_obj.ext_rem_im), 5, bub_im)
        
        bub_pixels = np.sum(bub_im == 10)
        gal_pixels = np.sum(bub_im == 0)
        
        # finding actual bubble coverage of galactic disk
        bub_gal_fraction = np.round(bub_pixels/(bub_pixels+gal_pixels), 3)

        # creating proper physical size axes tick marks
        edge_kpc = alg_obj.im_size/2
        nticks = 9
        xlabel = np.round(np.linspace(-edge_kpc, edge_kpc, nticks), 1)
        ylabel = np.round(np.linspace(edge_kpc, -edge_kpc, nticks), 1)
        ticks = np.round(np.linspace(0, alg_obj.im.shape[0], nticks), 0)

        # plotting bubbles or regions with custom colormap:
        #   - purple: galaxy region (not a bubble)
        #   - maroon: exterior region
        #   - yellow: detected bubbles or full regions (not masked based on size)
        fig, ax = plt.subplots(1, 1)
        cust_cmap = ListedColormap(["purple", "maroon", "yellow"])
        ax.imshow(bub_im, cmap = cust_cmap)

        # manipulating plot tick marks and labels
        ax.set_xticks(ticks)
        ax.set_yticks(ticks)
        ax.tick_params(axis='both', which='both', colors='white')
        ax.set_xticklabels(xlabel, rotation = 33, color = "black")
        ax.set_yticklabels(ylabel, color = "black")
        
        # cleaning up final plot
        plt.tight_layout()
        
        # returning bubble image forexterior labeling and saving
        return ax, bub_gal_fraction