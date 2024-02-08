'''File contains all of the functions required to apply Johnson-Nyquist Noise to a simulation 
galaxy given in H/cm^2 column density array form, as well as the functions needed to apply
a beam gaussian to a given image.'''

# imports
import numpy as np
from scipy.fft import fft2, ifft2, fftfreq
from scipy.stats import norm
import scipy.stats as stats


def take_fft(sim, obs, real = True):
    '''Function to take fast fourier transform of an image.
    
    Variables:
        - sim: image of the simulation that needs to have noise applied
        - obs: image of the observation whose noise will be used to apply
               noise to the simulation
        - real: boolean to choose whether the real or imaginary part
                of the fft is returned
    
    Output:
        - sim_fft: real or imaginary fft of the simulation image
        - obs_fft: real or imaginary fft of the observation image
        '''
    
    # checking if image is square
    if(np.shape(sim)[0] != np.shape(sim)[1]):
        print("Input Image must be Square")
        return
    
    # checking that images are the same size
    if(np.shape(sim) != np.shape(obs)):
        print("The images must be the same shape")
        return
    
    # taking fft of simulation and observation
    sim_fft = fft2(sim)
    obs_fft = fft2(obs)
    
    # splitting real and imaginary components based on condition
    if(real == True):
        sim_fft = np.real(sim_fft)
        obs_fft = np.real(obs_fft)

    elif(real == False):
        sim_fft = np.imag(sim_fft)
        obs_fft = np.imag(obs_fft)
    
    return sim_fft, obs_fft

def find_noise(obs_fft, seed = 123456789):
    '''Function to find the Johnson-Nyquist noise associated with a given 
    observation image.
    
    Variables:
        - obs_fft: 
        
    Output:
        - 
        '''
    
    # flatten to 1D array
    flat = obs_fft.flatten()
    N = np.shape(obs_fft)[0]
    
    # finding standard deviation of noise
    mean = np.mean(flat)
    x_adj = flat - mean
    x_sum = 0
    
    for i in range(len(x_adj)):
        x_sum += x_adj[i]**2
    std = np.sqrt(x_sum/(len(x_adj)-1))

    # calculating noise gaussian centered at 0 with given seed
    np.random.seed(seed)
    noise_gauss = np.random.normal(loc = 0, scale = std, size = (N, N))

    return noise_gauss

def recombine_parts(real_fft, imag_fft):
    # finding shape of image
    N = np.shape(real_fft)[0]
    
    # initializing output array
    arr_fft = np.zeros((N, N), complex)

    # recobining real and imaginary parts of fft array
    arr_fft.real = real_fft
    arr_fft.imag = imag_fft
    
    return arr_fft

def frequency_condition(arr_fft):
    # finding shape of image
    N = np.shape(arr_fft)[0]
    
    # finding the frequencies as determined by the FFT
    kfreq = fftfreq(N)*N
    kfreq2D = np.meshgrid(kfreq, kfreq)
    comb_freq = np.zeros((N,N,2))

    # organizing the kx and ky frequencies into one array paired at each pixel postion
    for i in range(len(kfreq2D[0])):
        for j in range(len(kfreq2D[1])):
            kx = kfreq2D[0][i,j]
            ky = kfreq2D[1][i,j]
            comb_freq[i,j] = np.array([kx, ky])
    
    # determining min and max frequency values
    min_kx = int(np.min(kfreq2D[0]))
    max_kx = int(np.max(kfreq2D[0]))
    min_ky = int(np.min(kfreq2D[1]))
    max_ky = int(np.max(kfreq2D[1]))
    
    # checking I(w) = I(-w)* condition for noise fft array
    for i in range(min_kx, max_kx + 1):
        for j in range(0, max_ky + 1):
            # finding corresponding fft array values
            pos_val = arr_fft[j, i]
            n_val = arr_fft[-j, -i]
            
            # forcing the frequency if not met
            if(pos_val != np.conj(n_val)):
                arr_fft[-j, -i] = np.conj(pos_val)
    
    # accounting for most 0 and nyquist frequencies (four edges minues the corners)
    for p in [0, max_kx+1]:
        for k in range(-max_kx, max_kx+2):
            arr_fft[p, k] = np.conj(arr_fft[p, -k])
            arr_fft[k, p] = np.conj(arr_fft[-k, p])

    # forcing the corners to be the absolute value since a real number is conjugate a real number                
    arr_fft[0,0] = np.abs(arr_fft[0,0])
    arr_fft[max_kx+1,max_kx+1] = np.abs(arr_fft[max_kx+1,max_kx+1])
    arr_fft[max_kx+1,0] = np.abs(arr_fft[max_kx+1,0])
    arr_fft[0,max_kx+1] = np.abs(arr_fft[0,max_kx+1])
    
    # checking all elements to see if frequency condition holds
    for i in range(min_kx, max_kx+1):
        for j in range(min_kx, max_kx+1):
            if(arr_fft[j,i] != np.conj(arr_fft[-j,-i])):
                print("\nProblem at\ni=", i, "\nj=", j)
    
    return arr_fft

def sim_noise_alg(sim, obs, seed = 123456789):
    # finding fft arrays
    real_sim, real_obs = take_fft(sim, obs, real = True)
    imag_sim, imag_obs = take_fft(sim, obs, real = False)
    
    # recombining sim components
    sim_fft = recombine_parts(real_sim, imag_sim)
    
    # finding noise
    real_noise = find_noise(real_obs, seed)
    imag_noise = find_noise(imag_obs, seed)
    
    # recombining noise components
    noise_fft = recombine_parts(real_noise, imag_noise)
    
    # correcting noise frequency to condition
    corr_noise_fft = frequency_condition(noise_fft)
    
    # adding noise to simulation
    sim_noise_fft = sim_fft + corr_noise_fft
    sim_noise = ifft2(sim_noise_fft)
    
    # checking imaginary component and removing if not large
    if(np.max(np.imag(sim_noise)) > 1e8):
        print("Imaginary components are too large and most likely were caused by an error in input.")
        return
    sim_noise = np.real(sim_noise)
    
    return sim_noise, sim_noise_fft

def fix_im(im):
    '''Correct Image after applying a Smoothing Gaussian.
    
    Variables:
        - im: the input image we are trying to correct
        
    Retruns:
        - out_im: the corrected image
    '''
    
    # finding size of image
    N_beam = im.shape[0]
    fix_N = int(N_beam/2)

    # seperating image into quadrants to rearrange
    tl = im[fix_N:, fix_N:]
    tr = im[fix_N:, :fix_N]
    bl = im[:fix_N, fix_N:]
    br = im[:fix_N, :fix_N]

    # recombining the image
    top_im = np.concatenate([tl, tr], axis = 1)
    bottom_im = np.concatenate([bl, br], axis = 1)
    out_im = np.concatenate([top_im, bottom_im], axis = 0)
    
    return out_im

def apply_beam(axes, N, pixel_size, im):
    '''Find and Apply Gaussian of the Beam Size.
    
    Variables:
        - axes: list of major and minor axes of beam in arcseconds [major, minor]
        - N: number of pixels on one side of image, must be square
        - pixel_size: the angular size of each pixel in arcseconds
        - im: the image that the beam will be applied to
        
    Returns: 
        - im_beam: the given image with the beam applied
     '''

    # beam dimensions
    if(len(axes) != 2):
        print("The beam axes must be given in a 2 sized list format.")
    major = np.max(axes)   # arcsecs
    minor = np.min(axes)   # arcsecs

    # image length in arcseconds
    im_len = pixel_size * N
    x_im = np.linspace(-0.5*im_len, 0.5*im_len, N)
    y_im = x_im
    X_im, Y_im = np.meshgrid(x_im, y_im)

    # calculating the beam std and gaussian
    beam_std_major = major/np.sqrt(8*np.log(2))
    beam_std_minor = minor/np.sqrt(8*np.log(2))
    beam_gauss = (1/(2*np.pi*beam_std_major*beam_std_minor))*np.exp(-0.5*((X_im/beam_std_major)**2 +(Y_im/beam_std_minor)**2))

    # applying beam to image in fourier space
    im_fft = fft2(im)
    beam_fft = fft2(beam_gauss)
    im_beam_fft = im_fft*beam_fft
    im_beam = ifft2(im_beam_fft)

    # checking that imaginary component is not large then removing imaginary components
    if(np.max(np.imag(im_beam)) > 1e8):
        print("Imaginary components are too large and most likely were caused by an error in input.")
        return
    im_beam = np.real(im_beam)

    # correcting the smoothed image
    im_beam = fix_im(im_beam)
    
    im_beam = np.where(im_beam < 1e18, 1e18, im_beam)
    im_beam = np.where(im_beam > 1e23, 1e18, im_beam)

    return im_beam