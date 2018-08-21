from StarCluster import *

'''
 This is an example of how to draw an image and search for sources (star candidates)
 You can run this for your images several times with different parameters
 until you optimize the star dentification, the cluster is centered in the area and
 you have identified one or more reference stars

'''

areaCluster = [800, 1200, 800, 1200]
fr = StarCluster('data/M52_example/', area = areaCluster)
fr.RemoveBkgImages()
fr.SetBkgSigmaThr(4)
fr.SetSourceFWHM(5)
fr.SetSourceRadius(10)

area = [0, 3000, 0, 3000]
fr.GetListOfBrightestSources('B', 20, area)
fr.DrawImage('B', '', 'sources', 20, area)

