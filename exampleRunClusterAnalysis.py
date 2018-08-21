from StarCluster import *
from StarPlotter import *

# Set path to FITs and output directoty, name of the cluster
path = 'data/M52_example/'
outdir = './output_example_M52/'
name = 'M52'

# Select area of the image with the object
clusterArea  = [1280, 1710, 960, 1380]

# Set distance of the cluster
distInParsec = 153. 

# Set reference stars
refStars = []
refStars.append(Star([157 , 1035], name = 'HIP 115198', R = -99, V = 6.944, B = 6.948, distance = lyToParsec(1664) ) )
refStars.append(Star([1464, 1531], name = 'HIP 115661', R = -99, V = 7.85, B = 8.56))
refStars.append(Star([1431,   46], name = 'HIP 115218', R = -99, V = 6.40, B = 8.07, distance = lyToParsec(1058) ))
refStars.append(Star([1398, 1087], name = 'HIP 115542', R = 7.84, V = 8.292, B = 9.345))

### Create the StarCluster area
fr = StarCluster(path, area = clusterArea, refStar = refStars, outdir = outdir, clusterName = name, verbose = 1)

# Automatically remove the background
fr.RemoveBkgImages()

# Tell if you want a RGB to BVR transformation
fr.SetRGBtoBVR(True)

# Set some parameres
fr.SetBkgSigmaThr(4)
fr.SetSourceFWHM(5)
fr.SetSourceRadius(10)
fr.SetDistance(distInParsec)

# Run the star identification and calibration
fr.SetStars()
fr.CalibrateStars()

# Save all the info
fr.DrawAll()

# Create a plotter an save some diagrams
plotter = StarPlotter(name, './M52/', addTaxis = True, baseStarArea = 150)
plotter.CreateColorMagDiagram('ColorMagnitude')
plotter.CreateHRdiagram('HRdiagram')
plotter.DrawRmVtoBmV('RmVtoBmV')
