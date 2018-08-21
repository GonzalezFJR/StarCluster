# StarCluster

### What is this?
  
  This repository contains some simple tools to identify stars form FITs images, calibrate the magnitudes and obtain some physical properties. It also prints data tables and produce some plots. 

  The code is based on astropy and photoutils.

  Requirements:
  - python 2.7 (do not know if it works in python 3.X)
  - numpy
  - matplotlib
  - astropy
  - photoutils

### Download the code

  git clone https://github.com/GonzalezFJR/StarCluster.git

### Images and other inputs

  You need images in FITs format. The header must contain the information of each filter.
  The filters must be: 'L', 'R', 'G' or 'V', 'B', and the last three are mandatory (RGB).
  The images must be aligned. This code do not align or stack images (yet).

  You must set at least one reference star, with magnitudes for each filter. You have to
  provide the position (in pixels) of that star in your images.

### Sumary on contents:

  - Several functions to load images, get sources, etc
  - Several physic functions and parametrizations
  - Class 'Star'
  - Class 'StarCluster'
  - Class 'StarPlotter'
  - Example of use (with some FITs files)

### Run the example

  First, run the script exampleDrawImg.py. This will draw an image for a given filter with 
  the identified sources. The 20 brightest sources are shown and their coordinates are printed.

  You can zoom in and search for some reference stars between the brightest stars.
  You can select the area of search of brightest stars or the area where the cluster should be.

     python exampleDrawImg.py

  After, you can run the example analysis. The output should be in a folder 'output_example_M52/'

     python exampleRunClusterAnalysis.py

### Todo list

  - Apply reddening corrections
  - Size of the points proportional to star radii
  - Better color template for temperatures
  - Add functions to StarPlotter with interesting phisics plots
  - Add plots with calibration for ref Stars
