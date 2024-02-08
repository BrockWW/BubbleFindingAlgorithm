'''File to convert 21 cm line emission radio astronomy observations to column density in
H*cm^-2. Input file to do conversion must be FITS file in units of Jy * beam^-1 * m * s^-1, 
example of which are images from THINGS database.'''

# imports
import numpy as np
from astropy.io import fits
from yt_to_column_density_functions import *

################################################################
# required values for all functions to run
################################################################

# file path of FITS file for conversion
obs_file_path = "/user/example/filename.FITS"

# file path to save output column density array to (without file extention)
save_file_path = "/user/example/filename"

################################################################
# creating conersion algorithm from equations (1) and (5) from 
# Walter et al.
################################################################

# obs data column density
obs_data = fits.getdata(obs_file_path, ext=0)[0][0]   # Jy * beam^-1 * m * s^-1   # from header of FITS file

print("Shape of observation image is ", obs_data.shape)

# conversion for NGC6946 equations (1) and (5) from
# Walter, F., Brinks, E., de Blok, W. J., Bigiel, F., Kennicutt, R. C., 
#     Thornley, M. D., & Leroy, A. (2008). Things: The H I nearby galaxy survey. The Astronomical Journal, 136(6), 2563–2647.
#     https://doi.org/10.1088/0004-6256/136/6/2563 

# eq (1)
TB_delV = (6.07e5) * obs_data / (6.04 * 5.61) / (1e3)   # K*km*s^-1
# eq (5)
obs_colden = (1.823e18) * TB_delV   # H*cm^-2
del(obs_data, TB_delV)

# setting minimum observation image value to remove small and negative values
obs_colden = np.where(obs_colden <= 1e18, 1e18, obs_colden)

# saving file for use in algorithm later
save_array(obs_colden, save_file_path)

print("Script complete.")