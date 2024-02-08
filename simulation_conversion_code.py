'''File combining all functons in order to convert hdf5 particle simulation
to atomic and molecular hydrogen column density arrays. Alter variables in
first block of code to customize output.'''

# imports
import sys
import numpy as np
from yt_to_column_density_functions import *

################################################################
# required values for all functions to run
################################################################

# z, y, z dimensions of 3D arrays from simulation
arr_d = [1024, 1024, 1024]
# physical size of x, y, z of simulation (pc)
dta = [40000, 40000, 20000]

# file path and name of hdf5 simulation file
file_p = "/user/example/filename.hdf5"

# file path and name (no file extension) to save density and 
# temperature 3D arrays
den_path = "/user/example/filename"
temp_path = "/user/example/filename"

# file path and name (no file extension) to save atomic and molecular 
# hydrogen column density arrays
atomicH_path = "/user/example/filename"
H2_path = "/user/example/filename"

################################################################
# checking memory before continuing script
################################################################

# calculating memory that saved
memory = (16/2097152) * np.prod(arr_d)
print("The given dimension parameters will produce arrays that require", memory, "MB")
print("The total memory that is required to successfully run the program is", 7*memory, "MB")

c = input("With the required memory, do you want to continue? (y/n) ")

# checking to continue the program or stop
if(c == "y"):
    pass
else:
    sys.exit("Program successfully aborted.")

################################################################
# calculating density and temperature from hdf5 simulation file
################################################################

den, temp = temp_den_arr(file_p, arr_d, dta)

print("Density and temperature converted to arrays.")

save_array(den, den_path)
save_array(temp, temp_path)

del(den, temp)

################################################################
# converting density to atomic and molecular hydrogen then 
# column density and saving
################################################################

density_atomicH, density_H2 = calculate_atomicH(den_path, temp_path, arr_d, dta)

print("Density converted to atomic and molecular hydrogen.")

colden_atomicH = colden_conversion(density_atomicH, N = arr_d[2], dl = dta[2], axis_num = 2)
colden_H2 = colden_conversion(density_H2, N = arr_d[2], dl = dta[2], axis_num = 2)

print("Conversion to column density for both atomic and molecular hydrogen complete.")

del(density_atomicH, density_H2)

save_array(colden_atomicH, atomicH_path)
save_array(colden_H2, H2_path)

print("Script complete.")