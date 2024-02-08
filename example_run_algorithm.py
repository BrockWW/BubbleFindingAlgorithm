'''A'''

# imports
from bubble_algorithm_class import *
from bubble_plotting_functions import *
import matplotlib.pyplot as plt
import numpy as np



################################################################################################
# Variables to change based on galaxy images used
################################################################################################
'''Variables:
    - seed: random seed for numpy to keep noise algorithm results the same (optional)
    - smooth_ims: boolean to decide whether to apply a smoothing gaussian to images based
                  on the observation galaxy's beam size
    - smooth: multiple of observation image's beam used to create smoothing gaussian
    - pixel_size: angular pixel side length for observation galaxy image (arcsec)
    - beam_size: beam semimajor and semiminor axis from observation galaxy image (arcsec)
    - obs_gal_dist: physical distance to observation galaxy (kpc)
    - sim_diam_kpc: physical size of simulationgalaxy image if used (kpc)
    - obs_file: file path of saved observation column density data .npy file (exlude extension)
    - obs_file: file path of saved simulation column density data .npy file (exlude extension)
    '''

seed = 92947490
smooth_ims = True
smooth = 2
pixel_size = 1.5   # arecseconds
beam_size = np.array([6.04, 6.04])   # arcseconds
obs_gal_dist = 5.9*1000   # kpc
sim_diam_kpc = 40   # kpc

obs_file = "/user/example/filename"
sim_file = "/user/example/filename"

################################################################################################
# Loading in galaxy images
################################################################################################
obs_im = np.load(obs_file + ".npy")   # H*cm^-2
sim_im = np.load(sim_file + ".npy")   # H*cm^-2

################################################################################################
# Finding physical size of each real galaxy
################################################################################################

# angular to physical size formula for some real galaxy
obs_pixel_size = pixel_size / 206265  # arcsec -> rad
obs_diam_kpc = (obs_pixel_size*obs_im.shape[0]) * obs_gal_dist   # kpc

################################################################################################
# Applying the observation image's beam and noise to the simulation image
################################################################################################

sim_beam = apply_beam(beam_size, sim_im.shape[0], pixel_size, sim_im)   # sim image with beam
sim_noise_beam, _ = sim_noise_alg(sim_beam, obs_im, seed)   # sim image with noise and beam

################################################################################################
# Using the Apply Beam algorithm to smooth simulation and observation images (optional)
################################################################################################

if(smooth_ims == True):
    # smoothing each galaxy with a gaussian of size a multiple of observation's beam size
    smoothing_axes = smooth*beam_size
    out_obs_im = apply_beam(smoothing_axes, obs_im.shape[0], pixel_size, obs_im)
    out_sim_im = apply_beam(smoothing_axes, sim_noise_beam.shape[0], pixel_size, sim_noise_beam)

else:
    out_obs_im = obs_im
    out_sim_im = sim_noise_beam

################################################################################################
# Running bubble finding algorithm for both galaxies and the 50-80th volume percentile
# and plotting resulting bubbles
################################################################################################

# selecting range of percentiles to analyze
step = 5
percentile_arr = np.arange(50, 85, step)

# collecting all needed data into lists for easy looping
name_list = ["obs_name", "sim_name"]
gal_ims = [out_obs_im, out_sim_im]
diam_list = [obs_diam_kpc, sim_diam_kpc]

# selecting correct variables
for i in range(len(gal_ims)):
    name = name_list[i]
    gal_use = gal_ims[i]
    gal_diam = diam_list[i]
    
    # selecting volume percentile to use
    for pe in percentile_arr:
        galaxy_obj = Bubble_Alg(gal_use, gal_diam, pe, beam_size[0], pixel_size)
        
        # running bubble plotting funtion
        bubble_ax, bub_gal_frac = ext_bubble_plotting(galaxy_obj)
        log_thresh = np.log10(galaxy_obj.colden_thresh)

        # properly labeling image
        bubble_ax.set_xlabel("x $(kpc)$")
        bubble_ax.set_ylabel("y $(kpc)$")
        bubble_ax.set_title("{n} Bubble Image\nVolume Percentile {p}".format(n = name, p = pe))

        # including bubble covereage fractions and colume density threshold onto plot
        bubble_ax.annotate("Log10 Threshold:\n{l}".format(l = np.round(log_thresh, 3)), (0.6, 0.88), xycoords = "axes fraction", **{"color": "white", "fontsize": "small"})
        bubble_ax.annotate("Bubble-Galaxy\nPixel Ratio:\n{r}".format(r = bub_gal_frac), (0.1, 0.83), xycoords = "axes fraction", **{"color": "white", "fontsize": "small"})

        # saving and closing image
        plt.savefig("TEST_IMS/{n}_{p}.png".format(n = name, p = pe), bbox_inches = "tight")
        plt.clf()
        plt.close()