Title and Description
=====================

Superbubble Finding Algorithm

We create an algorithm that is able to automatically detect and measure superbubbles left in
galactc disks using common python libraries. The algorithm takes an input column density galaxy image
and returns the labels and basic measurements of the detected bubbles so tasks like plotting
and further analysis can be done. Functions to solve for the superbubble radii, superbubble galactic 
radii location, and an external function to plot the detected bubbles are prepackaged in the algorithm.

#############################################################################################################

Contained Files
===============

    - bubble_algorithm_class.py: Contains the main part of the algorithm where a class structiure is used
                                 to conviniently store and access different values during analysis.
    - bubble_plotting_functions.py: Contains a function used to plot the detected superbubbles after running 
                                    the algorithm.
    - example_run_algorithm.py: An example of how this algorithm should be used and what it requires to
                                run properly. Please refer to this file to understand general setup and use.
    - noise_algorithm: Contains required functions for applying noise and beam size to a simulation image
                       based on a radio observation image's noise and beam.
    - radio_observation_conversion.py: Contains the code needed to convert FITS radio astronomy images in
                                       units of Jy * beam^-1 * m * s^-1 to the required units of H*cm^-2.
                                       (Not required)
    - README.txt: This file.
    - required_libraries.txt: File containing the necessary Python libraries and their working versions to
                              run the algorithm. See installation and use notes below.
    - required_libraries_noversion.txt: File containing the necessary python libraries without versions to
                                        run the algorithm. See installation and use notes below.
    - simulation_conversion_code.py: File to convert particle .hdf5 galaxy simulation data files into readable
                                     numpy arrays in correct units for algorithm. (Not required)
    - test_data: Simulation and observation data included for testing the algorithm.
    - test_ims: Directory to store output images created from running example_run_algorithm.py using the files
                provided in test_data.
    - yt_to_column_density_functions.py: File containing the necessary functions in order to convert the units
                                         of the array simulation data to units of H*cm^-2. (Not required)

#############################################################################################################

Installation and Use
====================

To import this algorithm, clone the github repository using: 

git clone https://github.com/BrockWW/BubbleFindingAlgorithm.git

This algorithm has been developed and tested using Python 3.6.4, all the required libraries and versions needed
specifically for Python 3.6.4 are described in required_libraries.txt. If using a more recent version of python, 
please use the required_libraries_noversion.txt to correctly install all necessary python libraries. The algorithm 
has been tested with Python 3.12 using the most up to date libraries and gives results as was seen with Python 3.6.4. 
Please ensure that the local python installation being used is allowed to overwrite currently installed libraries. 
It is recommended a virtual environment is created with the libraries provided to prevent any file conflict. To 
install all of these libraries, using the terminal move into the location of the imported algorithm code and run the 
following command:

pip install -r required_libraries.txt

or

pip install -r required_libraries_noversion.txt

Run the first "pip install" line if Python 3.6.4 is being used or run the second line if a differnt Python version is used.
This will install all of the listed libraries and their correct versions to the local python installation. From here,
it is recommended that the file example_run_algorithm.py be used as a reference when running the algorithm as it goes
through the steps of reading in files, to running the algorithm and plotting the output. The example data file type used is
.npy, but as long as the input into the algorithm is a 2D numpy array in units of H*cm^-2, any type of file can 
be used for importing the data.

To run the files, simply use the normal Python routine in terminal of typing:

python FILE_NAME.py

**NOTE: The files:
    - example_run_algorithm.py
    - radio_observation_conversion.py
    - simulation_conversion_code.py
all require some user manipulation to run properly. Most of this includes specifying the path of locally stored files.
If no change is made to match the user's file locations and save directories, these will not run and the algorithm will fail.