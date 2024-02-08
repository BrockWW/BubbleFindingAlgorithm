'''File containing all the required functions to convert an HDF5 particle galaxy simulation into
column density arrays of atomic and molecular hydrogen. Because of large memory requirements and
long running time, it is highly recommended that arrays from steps be saved to .npz files so if
a later step fails, the entire process does not need to be redone.'''

# imports
import numpy as np
import yt

def save_array(arr, file_path):
    '''Function to save any given numpy array to a specific path location as an
    .npy file.
    
    Variables:
        - arr: any numpy array of any size
        - file_path: file path for the location that the array will be saved 
                     to as an .npy file
            
    Output:
        - saves array as .npy file
         '''
    
    # check if path is given as string
    if(isinstance(file_path, str)):
        np.save(file_path, arr)
        return
    
    # return error if path not given as string
    else:
        print("Path needs to be string with file name, discluding file extension.")
        return

def temp_den_arr(file_path, arr_dims = [1024, 1024, 1024], dtau = [40000, 40000, 20000]):
    '''Function to convert baryon particle HDF5 galaxy simulation density and temperature 
    values to 3D numpy arrays of a given size. It is highly recommended to have the same
    x and y dimensions, as later algorithm functions rely on the input images being square.
    
    Variables:
        - file_path: location of HDF5 file
        - arr_dims: list of number of cells in each [x,y,z] dimension
        - dtau: list of physical size of galaxy in each [x,y,z] dimension (pc)
           
    Output:
        - density: 3D array of density saved from HDF5 file (MSun*pc^-3)
        - temperature: 3D array of temperature saved from HDF5 file (K)
        '''
    
    # checking length of arr_dims and dtau
    if(len(arr_dims) != 3 or len(dtau) != 3):
        print("An array of length 3 must be given specifying the\
              [x,y,z] dimensions wanted for the output array for\
              both arr_dims and dtau.")
        return
    
    # set new unit base for yt to read in
    units_override = {"UnitLength_in_cm": 3.085678e21,
                      "UnitMass_in_g": 1.989e43,
                      "UnitVelocity_in_cm_per_s": 1e5}
    
    # check if path is given as string
    if(isinstance(file_path, str)):
        # load in the dataset
        ds = yt.load(file_path, unit_base = units_override)
        
        # finding half of size of galaxy in kpc
        dx = dtau[0]/2000
        dy = dtau[1]/2000
        dz = dtau[2]/2000

        # define an arbitrary grid to move particles types to cell types
        grid_obj = ds.arbitrary_grid([-dx, -dy, -dz], [dx, dy, dz], dims=arr_dims)

        # deposit the particles onto grids for the following variables
        density = np.array(grid_obj["PartType0", "Density"].in_units("Msun/pc**3"))
        temperature = np.array(grid_obj["PartType0", "Temperature"].in_units("K"))
        
        return density, temperature
    
    else:
        print("Input file path needs to be a string")
        return
        
def calculate_atomicH(density_file, temp_file, arr_dims = [1024, 1024, 1024], dtau = [40000, 40000, 20000]):
    '''Algorithm created using equations from:
    Krumholz, M. R., & Gnedin, N. Y. (2011). A comparison of methods for determining the molecular 
        content of Model Galaxies. The Astrophysical Journal, 729(1), 36. 
        https://doi.org/10.1088/0004-637x/729/1/36
        
    NOTE: Depending on the size of the density and temperature files, this function can use up a lot of memory
    as it creates new arrays before it is able to remove unneeded arrays. For whatever the size of the density
    or temprature array, this function will need at least 5 times the memory of one of the arrays.

    Variables:
        - density_file: path of file storing density array (MSun*pc^-3)
        - temp_file: path of file storing temperature array (K)
        - arr_dims: integer or list of number of cells in each x,y,z dimension
        - dtau: list of physical size of galaxy in each [x,y,z] dimension (pc)
    
    Returns:
        - density_atomicH: 3D array of the calculated density of atomic hydrogen (g*cm^-3)
        - density_H2:3D array of the calculated density of molecular hydrogen (g*cm^-3)
        '''
    
    # checking length of arr_dims
    if(len(arr_dims) != 3 or len(dtau) != 3):
        print("An array of length 3 must be given specifying the\
              [x,y,z] dimensions wanted for the output array for\
              both arr_dims and dtau.")
        return

    
    '''Initialization of important arrays and variables.'''
    # finding the length of a cell in the x, y, and z direction in cm
    dx = (dtau[0]/arr_dims[0]) * (3.086e18)   # cm
    dy = (dtau[1]/arr_dims[1]) * (3.086e18)   # cm
    dz = (dtau[2]/arr_dims[2]) * (3.086e18)   # cm

    # importing saved data from files
    density_use = np.load(density_file + ".npy")   # MSun/pc^3
    temp_use = np.load(temp_file + ".npy")   # K
    
    print("den shape", density_use.shape)
    print("temp shape", temp_use.shape)

    # checked this, correctly reshapes data
    #density_use = density_use.reshape(arr_dims[0], arr_dims[1], arr_dims[2])
    #temp_use = temp_use.reshape(arr_dims[0], arr_dims[1], arr_dims[2])

    # applying mass fraction of Hydrogen in galaxy
    XH = 0.71
    density_use *= XH

    # overwritting the zero values that appear by using small
    # values for temperature and density
    density_use[np.where(density_use == 0)] = 1e-6
    temp_use[np.where(temp_use == 0)] = 1

    
    '''Code to find gradient in units g/cm^4.'''
    density_noionized = np.where(temp_use < 1e4, density_use, 1e-6)   # MSun/pc^3
    del(density_use)
    del(temp_use)

    den_conv = (1.9891e33) / (2.93799895e55)   # g*pc^3/Msun/cm^3
    H_den = density_noionized * den_conv   # g/cm^3
    del(density_noionized)

    grad_rho_p = np.gradient(H_den, dx, dy, dz)   # g/cm^3/cm
    grad_rho = np.sqrt(grad_rho_p[0]**2 + grad_rho_p[1]**2 + grad_rho_p[2]**2)   # magnitude of gradient
    del(grad_rho_p)

    
    '''Code to find scale height in units cm.'''
    mu_H = 2.3e-24   # grams
    sigma_d = 1e-21   # cm^2

    h = H_den/np.abs(grad_rho)   # cm
    h = np.where(h > (100*3.086e18), (100*3.086e18), h)  # set max h value to 100pc
    del(grad_rho)

    
    '''Code to find unitless optical depth.'''
    eps = H_den * h   # g/cm^2
    del(h)
    
    tau_C = (eps*sigma_d)/mu_H   # unitless
    del(eps)

    
    '''Code to calculate fraction of H2 and atomic H.'''
    X = 71/(H_den/(1.6735575e-24))   # divide by number of H per cm^3, eq. (3)

    s = np.log(1 + (0.6*X) + (0.01*X**2))/(0.6*tau_C)   # eq. (2)
    del(X)
    del(tau_C)

    fH2 = 1 - ((3/4)*(s/(1 + 0.25*s)))   # approx function, eq. (1)
    del(s)
    fH2 = np.where(fH2 < 0, 1e-10, fH2)   # Forcing minimum to be ~0
    
    fH = 1 - fH2

    density_atomicH = H_den * fH   # g/cm^3
    density_H2 = H_den * fH2   # g/cm^3
    del(fH)
    del(fH2)
    del(H_den)
    
    return density_atomicH, density_H2

def colden_conversion(den_arr, N = 1024, dl = 20000, axis_num = 2):
    '''Function that converts a 3D density array given in g*cm^-3 into a 2D column 
    density array in H*cm^-2.
       
    Variables:
        - den_arr: 3D density array to convert to column density (g*cm^-3)
        - N: number of cells in the axis that the density will be summed in
        - dl: physical size of image in the axis that the density will be summed in (pc)
        - axis_num: the axis that will be used to sum the density over
       
    Outputs:
        - colden_arr: 2D array that gives the column density of the galaxy (H*cm^-2)
       '''
    
    # calculating the column density
    dl = (dl/N) * (3.086e18)   # cm
    colden_arr = np.sum(den_arr, axis = axis_num)*dl/(1.6735575e-24)   # H*cm^-2
    
    # relic of orienting simulation galaxy in numpy to same orientation from yt
    colden_arr = colden_arr.transpose()
    colden_arr = np.flip(colden_arr, axis=0)
    
    return colden_arr