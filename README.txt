Title and Description
=====================

Supernovae Bubble Finding Algorithm

We create an algorithm that is able to automatically detect and measure supernovae bubbles left in
galactc disks using common python libraries. The algorithm takes an input column density galaxy image
in units of hydrogen*cm^-2 and returns the labels of the found bubbles so other manipulations like plotting
or other external measurements can be completed. The algorithm comes with functions to solve for the bubble
radii and galactic radii, and an external function to plot the found bubbles in a separate image.

#############################################################################################################

Contained Files
===============

    - bubble_algorithm_class.py: Contains the main part of the algorithm where a class structiure is used
                                 to conviniently store and access different values during analysis.
    - bubble_plotting_functions.py: Contains a function used to plot the found bubbles after running the
                                    algorithm.
    - example_run_algorithm.py: An example of how this algorithm should be used, and what it requires to
                                run properly is included, refer to this to understand general seup and use.
    - noise_algorithm: Contains required functions for applying noise and beam size to a simulation image
                       based on an observation images noise and beam size.
    - radio_observation_conversion.py: Contains the code needed to convert FITS radio astronomy images in
                                       units of Jy * beam^-1 * m * s^-1 to the required units of H*cm^-2.
                                       (Not required)
    - README.txt: This file.
    - required_libraries.txt: File containing the necessary python libraries and their working versions to
                              run the algorithm. See installation and use notes below.
    - simulation_conversion_code.py: File to convert particle .hdf5 galaxy simulation data files into readable
                                     numpy arrays in correct units for algorithm. (Not required)
    - test_data: Simulation and observation data included for testing the algorithm. (To be removed)
    - yt_to_column_density_functions.py: File containing the necessary functions in order to do the simulation
                                         to H*cm^-2 array conversion. (Not required)

#############################################################################################################

Installation and Use
====================

To import this algorithm, clone the github repository using: 

git clone https://github.com/BrockWW/BubbleFindingAlgorithm.git

This algorithm has been developed and tested using Python 3.6.4, all the required libraries and versions needed 
are described in required_libraries.txt. Please ensure that the local python installation being used is allowed 
to overwrite currently installed libraries. It is recommended to create a virtual environment if it is not feasable 
to install this package with the local main python installation. In order to install all of these libraries, using 
the terminal move into the location of the imported algorithm code and run the following command:

pip install -r required_libraries.txt

or

pip install -r required_libraries_noversion.txt

This will install all of the listed libraries and their correct versions to the local python installation. From here,
it is recommended that the file example_run_algorithm.py be used as a reference when running the algorithm as it goes
through the steps of reading in files to running the algorithm. All of the examples read in files of the .npy type, but
as long as the input into the algorithm is a 2D numpy array in units of H*cm^-2, any type of file can be used for
importing the data.

To run the files, simply use the normal python routine in terminal of typing:

python FILE_NAME.py

**NOTE: The files:
    - example_run_algorithm.py
    - radio_observation_conversion.py
    - simulation_conversion_code.py
all require some user manipulation to run properly. Most of this included specifying the path of locally installed files.
If no change is made to match the user's file locations and save directories, these will not run and the algorithm will fail.





