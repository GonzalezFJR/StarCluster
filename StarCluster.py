'''
  @ autor: Xuan
  @ mail: gonzalezrodrigo@uniovi.es

  Conjunto de funciones para leer imagenes .FIT y obtener cantidades fotometricas
  Lee los .FIT, obtiene estrellas, permite obtener calibraciones, etc.
  Las imagenes deben estar:
    - Alineadas
    - Contener el filtro en la cabecera (L, R, G/V, B) (recomendable, aunque se pueden introducir a mano)
    - Para obtener todas las cantidades, deben existir al menos medidas para filtos R,G,B
    - Debe existir una estrella de referencia y se deben conocer sus datos y posicion en la imagen

  ========================================================================================================
  ========= Resumen de contenidos

  Conjunto de funciones para derivar cantidades fisicas y realizar parametrizaciones.
  Temperatura, luminosidad, magnitud absoluta, radio estelar, etcetera.

  Clase Star
    Input: magnitud para cada filtro (ya calibrada), posicion en la imagen...
    Metodos para calcular todas las cantidades fisicas
  
  ClaseStarCluster 
    Inputs:
      - Path a las imagenes (y filtros en caso de que no esten en las cabezeras)
      - Area de interes
      - Estrella de referencia (objeto Star), debe contener posicion y magnitudes
      - Distancia al cumulo
      - Nombre del cumulo
      - Otros parametros de analisis (fwhm de ajuste de las estrellas, umbral respecto el fondo, etc)

    Metodos:
      - Resta el fondo
      - Calibra magnitudes
      - Identifica estrellas y dibuja area de referencia, 'sources' identificadas, etc
      - Crea los objetos Star a partir de los diferentes filtros, calculando todas las magnitudes fisicas posibles
      - Produce imagenes y tablas de la identificacion y calibracion y tablas de datos como output del analisis

  Clase StarPlotter
    Inputs: path a la tabla de datos (output de StarCluster) y parametros de estilo
    Dibuja el diagrama HR y otras representaciones graficas

  ==============================================================================================

  REQUISITOS (software)
  - python 2.7 (no he probado si funciona en python 3.x)
  - numpy
  - matplotlib
  - astropy y sus requisitos
  - photoutils y sus requisitos

  ==============================================================================================

  >> Ejemplo para trabajar con 'sources' (estellas candidatas en un filtro). Puedes hacer esto
  >> para identificar el area del cumulo y la posicion (en pixeles) de las estrellas de referencia:

    areaCluster = [800, 1200, 800, 1200] # De referencia, iremosmo dificandolo
    fr = StarCluster('data/M29/', area = areaCluster)
    fr.RemoveBkgImages()
    fr.SetBkgSigmaThr(4)
    fr.SetSourceFWHM(5)
    fr.SetSourceRadius(10)
    fr.GetListOfBrightestSources('V', 20, [0,3000,0,3000])
    fr.DrawImage('V', '', 'sources', 20, [0,3000,0,3000])

  >> Ejemplo de analisis de un cluster:

     # Define los inputs, incluyendo al menos una estrella de referencia (con su posicion el la imagen):
     clusterName = 'M52'
     clusterArea = [1280, 1710, 960, 1380]
     path = 'data/M52/'
     path = './M52/'
     refStar = Star([1398, 1087], name = 'HIP 115542', R = 7.84, V = 8.292, B = 9.345)

     # Crea el objeto StarCluster con los inputs y establece los parametros segun tus necesidades:
     fr = StarCluster(path, area = clusterArea, refStar = refStars, outdir = outdir, clusterName = name, verbose = 1)
     fr.RemoveBkgImages()
     fr.SetRGBtoBVR(True)
     fr.SetBkgSigmaThr(4)
     fr.SetSourceFWHM(5)
     fr.SetSourceRadius(10)
     fr.SetDistance(153.) # Parsec

     # Ejecuta la identificacion de estrellas y lacalibracion
     fr.SetStars()
     fr.CalibrateStars()

     # Guarda todos los resultados
     fr.DrawAll()

  >> Dibuja bonitos diagramas a partir de las tablas de datos, output del analisis:

     plotter = StarPlotter(name, './M52/', addTaxis = True, baseStarArea = 150)
     plotter.CreateColorMagDiagram('ColorMagnitude')
     plotter.CreateHRdiagram('HRdiagram')
     plotter.DrawRmVtoBmV('RmVtoBmV')

'''

import numpy as np
import os, sys
from photutils import datasets
from photutils import DAOStarFinder
import matplotlib.pyplot as plt
import astropy
from matplotlib.pyplot import figure
import matplotlib.patches as patches
import matplotlib.pyplot as plt
from photutils import aperture_photometry, CircularAperture
from astropy.visualization import *

#=========== Read fits and extract info
#========================================================================

def LoadFitImage(path, area = []):
  ''' Opens a .fit image and returns the float 2D array and the filter (read from the header) '''
  if not os.path.isfile(path): 
    print 'ERROR: not found file ' + path
    return
  hdu = astropy.io.fits.open(path)[0]  #datasets.load_star_image()
  filter = ''
  if 'FILTER' in hdu.header: filter = hdu.header['FILTER']
  image = hdu.data.astype(float)
  if not len(area) > 0:
    return image, filter
  x0, x1, y0, y1 = area
  recImg = hdu.data[x0:x1, y0:y1].astype(float)
  return image, filter, recImg

def GetBkgMedianAndSigma(image):
  ''' Obtains the background as median of the image and the sigma '''
  from astropy.stats import mad_std
  bkg = np.median(image)
  bkg_sigma = mad_std(image)
  return [bkg, bkg_sigma]

def GetSources(image, fwhm = 5., thr = 3, dataformat = '%.2g'):
  ''' Find stars with a given fwhm and threshold '''
  bkg,bkgsigma = GetBkgMedianAndSigma(image)
  daofind = DAOStarFinder(fwhm = fwhm, threshold = thr*bkgsigma)
  sources = daofind(image)
  for col in sources.colnames: sources[col].info.format = dataformat
  return sources

def GetBrightestSource(a):
  ''' Find the brightest star in an image or a group of stars '''
  if   isinstance(a, astropy.table.table.Table): # Table of stars
    maxflux = max(a['flux'])
    minmag  = min(a['mag'])
    bs = a[0]
    for s in a:
      #if s['flux'] == maxflux: bs = s
      if s['mag']  == minmag: bs = s
    return bs
  elif isinstance(a, numpy.ndarray): # Image
    return GetBrightestSources(GetSources(a))
  else:
    print '[GetBrightestSources] ERROR: unknown type ', type(a)

def GetListOfBightestSources(sources, nSources = 10, area = []):
  x = sources['xcentroid']; y = sources['ycentroid']; coor = np.transpose([x,y])
  magb = sources['mag']
  c = [z for _, z in sorted(zip(magb,coor))]
  xb,yb = np.transpose(c[:])#[:nSources])
  magb = sorted(magb)[:]#[:nSources]
  x = []; y = []; mag = []
  if len(area) >= 4:
    x0, x1, y0, y1 = area
    for i in range(len(xb)):
      xs = xb[i]; ys = yb[i];
      if xs < x1 and xs > x0 and ys < y1 and ys > y0:
        x.append(xs); y.append(ys); mag.append(magb[i])
    x = x[:nSources]; y = y[:nSources]; mag = mag[:nSources]
  else:
    x = xb[:nSources]; y = yb[:nSources]; mag = magb[:nSources]
  return [x, y, mag]

def GetSourceInPos(sources, px, py, thr = 5):
  ''' Looks for a source in [px, py] from all stars in the given list '''
  for s in sources:
    x = s['xcentroid']
    y = s['ycentroid']
    if abs(x-px) < thr and abs(y-py) < thr: return s
  print '[GetSourceInPos] WARNING: not found source in position [x, y] = [%i, %i]'%(px, py)
  return 0

def GetStarInPos(stars, px, py, thr = 5):
  ''' Looks for a star in [px, py] from all stars in the given list '''
  mindist = 1000000
  for s in stars:
    dist = s.GetDistToPoint(px, py)
    if dist < mindist:
      mindist = dist
      closestStar = s
    if mindist < thr: return s
    else: print '[GetStarInPos] WARNING: not found star in position [x, y] = [%i, %i]'%(px, py)
  return 0

def GetMagnitudeCorrectionForStars(sources, xpos, ypos, mag):
  ''' Returns constant to calibrate mangitudes '''
  refStar = GetSourceInPos(sources, xpos, ypos)
  starMag = refStar['mag']
  k = mag + starMag
  return k

def GetMagCorrFromRefStar(sources, star, filter):
  ''' Returns constant to calibrate mangitudes using a reference star as input '''
  x = star.GetX()
  y = star.GetY()
  mag = star.GetMag(filter)
  return GetMagnitudeCorrectionForStars(sources, xpos, ypos, mag)

def CalibrateSources(sources, k):
  ''' Returns the list of sources with calibrated mangitude '''
  for s in sources: s['mag'] += k
  return sources

def CalibrateSourcesFromRefStar(sources, star, filter):
  ''' Returns the list of sources with calibrated magnitude using a reference star as input '''
  return CalibrateSources(sources, GetMagCorrFromRefStar(sources, star, filter))

def SetForm(val, n = 4, form = '%1.2f'):
  t = form%val
  while len(t) < n: t+=' '
  return t

def GetInfoStar(index, stars):
  ''' Gets all info from one star '''
  s = stars[index]
  t  = ''
  t += ' ' + SetForm(index,                 4, '%i')    + ' |'
  t += ' ' + SetForm(s.GetX(),              4, '%i')    + ' |'
  t += ' ' + SetForm(s.GetY(),              4, '%i')    + ' |'
  t += ' ' + SetForm(s.GetMag('B'),         6, '%1.3f') + ' |'
  t += ' ' + SetForm(s.GetMag('V'),         6, '%1.3f') + ' |'
  t += ' ' + SetForm(s.GetMag('R'),         6, '%1.3f') + ' |'
  t += ' ' + SetForm(s.GetAbsMag(),         6, '%1.3f') + ' |' 
  t += ' ' + SetForm(s.GetBolometricMag(),  6, '%1.3f') + ' |' 
  t += ' ' + SetForm(s.GetLumi(),           8, '%1.2g') + ' |' 
  t += ' ' + SetForm(s.GetBmV(),            6, '%1.3f') + ' |'
  t += ' ' + SetForm(s.GetT(),              5, '%1.0f') + ' |'
  t += ' ' + SetForm(s.GetRadius(),         5, '%1.2f') + ' |'
  t += '\n'
  return t

def PrintStarData(stars):
  t = ''
  head  = '\n'
  head += '======================================================================================================\n' 
  head += '  id  |  X   |  Y   |   B    |   V    |   R    | VMag   | BolMag | Lumi     |  B-V   | T (K) |  Rad   \n'
  head += '------------------------------------------------------------------------------------------------------\n' 
  #head+= ' 0001 | 9999 | 0000 | 00.000 | 00.000 | 00.000 | 00.000 | 00.000 | 1.00e+11 | 00.000 | 10000 | 10000  \n'
  for i in range(len(stars)):
    if i%50 == 0: t += head
    t += GetInfoStar(i, stars)
  return t



#=========== Formulas and parametrizations
#========================================================================

from numpy import log10

def GetBGRfromRGB(filter, R, G, B):
  ''' Gets the BGR transformation from RGB filters '''
  CB_BG = 0.280; CB_GR = 0.600;
  CV_BG = 0.542; CV_GR =-0.064;
  CR_BG = 0.051; CR_GR = 0.468;
  if   filter == 'B':
    f = B + CB_BG*(B-G) + CB_GR*(G-R)
  elif filter == 'V':
    f = B + CV_BG*(B-G) + CV_GR*(G-R)
  elif filter == 'R':
    f = B + CR_BG*(B-G) + CR_GR*(G-R)
  return f

def ApplyRGBtoBVRtoStar(star, globalfac = 0):
  oG = star.GetMag('V')+globalfac
  oR = star.GetMag('R')+globalfac
  oB = star.GetMag('B')+globalfac
  R = GetBGRfromRGB('R', oR, oG, oB)
  B = GetBGRfromRGB('B', oR, oG, oB)
  V = GetBGRfromRGB('V', oR, oG, oB)
  star.SetMag(R, 'R')
  star.SetMag(B, 'B')
  star.SetMag(V, 'V')

def Interpolate(T0, T1, k0, k1, T):
  ''' Interpolate between two values... used to get bolometric corrections '''  
  if float(T1)-float(T0) == 0.: return 0.
  f = (float(T) - float(T0))/(float(T1)-float(T0))
  k = k0 + (k1-k0)*f
  return k

def GetBolometricCorrection(T):
  ''' Returns a bolometric correction from filer V, based on the star temperature.
      Assumed a luminosity class V (mains sequence), but the correction is not so different for 
      Red giants or supergiants.
      Linear interpolation from tabulated data.
      Ref: http://www.as.utexas.edu/~sj/a358-sp06/lec3.pt2.pdf, Table 3.7. '''
  k = 0; T0 = 0; T1 = 0; k0 = 0; k1 = 1;
  if T < 2640:
    warningError = '[WARNING] No bolometric correction for Temperature < 2640 k!! That\'s too cold!!! (T = %1.0f k)';
    T1 =  2640; T0 =  0
    k1 = -4.10; k0 = -8
  elif T < 3240:
    T1 =  3240; T0 =  2640
    k1 = -2.73; k0 = -4.10
  elif T < 3850:
    T1 =  3850; T0 =  3240
    k1 = -1.28; k0 = -2.73
  elif T < 4350:
    T1 =  4380; T0 =  3850
    k1 = -0.72; k0 = -1.28
  elif T < 5250:
    T1 =  5250; T0 =  4350
    k1 = -0.31; k0 = -0.72
  elif T < 5770:
    T1 =  5770; T0 =  5250
    k1 = -0.21; k0 = -0.31
  elif T < 5860:
    T1 =  5860; T0 =  5770
    k1 = -0.20; k0 = -0.21
  elif T < 6030:
    T1 =  6030; T0 =  5860
    k1 = -0.18; k0 = -0.20
  elif T < 6440:
    T1 =  6440; T0 =  6030
    k1 = -0.14; k0 = -0.18
  elif T < 7200:
    T1 =  7200; T0 =  6440
    k1 = -0.09; k0 = -0.14
  elif T < 8200:
    T1 =  8200; T0 =  7200
    k1 = -0.15; k0 = -0.09
  elif T < 9520:
    T1 =  9520; T0 =  8200
    k1 = -0.30; k0 = -0.15
  elif T < 11900:
    T1 = 11900; T0 =  9520
    k1 = -0.80; k0 = -0.30
  elif T < 13000:
    T1 = 13000; T0 = 11900
    k1 = -1.02; k0 = -0.80
  elif T < 15400:
    T1 = 15400; T0 = 13000
    k1 = -1.46; k0 = -1.02
  elif T < 18700:
    T1 = 18700; T0 = 15400
    k1 = -1.94; k0 = -1.46
  elif T < 22000:
    T1 = 22000; T0 = 18700
    k1 = -2.35; k0 = -1.94
  elif T < 30000:
    T1 = 30000; T0 = 22000
    k1 = -3.16; k0 = -2.35
  elif T < 33000:
    T1 = 33000; T0 = 30000
    k1 = -3.33; k0 = -3.16
  elif T < 38000:
    T1 = 38000; T0 = 33000
    k1 = -3.68; k0 = -3.33
  elif T < 44500:
    T1 = 45000; T0 = 38000
    k1 = -4.40; k0 = -3.68
  else:
    T1 = 52500; T0 = 45000
    k1 = -4.75; k0 = -4.40
  k = Interpolate(T0, T1, k0, k1, T)
  return k

def GetTforBmV(d):
   ''' Parametric formula to get the effective temperature from B-V color index '''
   C1 = 3.979145;
   C2 = -0.654499;
   C3 = 1.74069;
   C4 = -4.608815;
   C5 = 6.7926;
   C6 = -5.39691;
   C7 = 2.19297;
   C8 = -0.359496;
   logt = C1 + C2*d + C3*d*d + C4 * d*d*d + C5 * d*d*d*d + C6 * d*d*d*d*d + C7 * d*d*d*d*d*d + C8 * d*d*d*d*d*d*d
   if(logt > 20): logt = 0
   return  pow(10,logt);

def plotRelationBmVvsT():
  ''' Plot the correction factor B-V to temperature  '''
  import matplotlib.pyplot as plt
  import numpy as np
  x = np.linspace(-0.3, 1.9, 1000)
  y = [GetTforBmV(ix) for ix in x]
  plt.plot(x,y,'r')
  plt.ylabel('Temperature (K)')
  plt.xlabel('B-V')
  plt.show()

def GetAbsMagnitudeFromMagAndDistance(m, d):
  ''' Returns the absolute magnitud from the relative magnitude and the distance '''
  # d in parsecs
  M = m - 5*(np.log10(d/10)-1)
  return M

def lyToParsec(d):
 ''' Lihgt years to parsecs '''
 return float(d)/3.2616

def GetBolometricMagnitude(M, T):
  ''' Returns the bolometric magnitud from the magnitude and the temperature '''
  Mbol = GetBolometricCorrection(T) + M
  return Mbol

def GetLumiFromMbol(mbol):
  ''' Gets the luminosity of a star in units of sun luminosity from the bolometric magnitude '''
  mbolSun = 4.74;
  l = pow(10, 0.4*(mbolSun - mbol))
  return l

def LumiSolToWatt(lum):
  ''' From solar luminisities to watts '''
  Lsun = 3.828e26 # 3.0128e28
  return lum*Lsun

def GetStarRadius(lumi, temp):
  ''' Gets the radius of the star using the total luminosity and Stefan-Botlzmann's law. 
      The input lumi must be in watt and temperature in kelvin.
      If lumi < 10^10, will assume it's on solar lumis '''
  if lumi < 1e10: lumi = LumiSolToWatt(lumi)
  sigma = 5.67e-8 # watt m-2 K-4
  T4 = temp*temp*temp*temp
  pi = 3.141592
  r2 = lumi/(4*pi*sigma*T4)
  solarRadius = 6.96e8
  return np.sqrt(r2)/solarRadius 


#===============================================================
#========== Drawing

def GetColor(t):
    ''' Get the color from a temperature '''
    if(t > 1000): t = float(t)/1000
    # < 2k: Brown
    if   t < 2: return '#630202'
    # 2k to 5.5k Red-orange-yellow
    elif t < 2.2: return '#d60000'
    elif t < 2.4: return '#f70c0c'
    elif t < 2.6: return '#fc351b'
    elif t < 2.8: return '#fc531a'
    elif t < 3.0: return '#ff6021'
    elif t < 3.2: return '#ff7220'
    elif t < 3.4: return '#ff8020'
    elif t < 3.7: return '#ff9620'
    elif t < 4.0: return '#ffb420'
    elif t < 4.3: return '#ffe120'
    elif t < 4.6: return '#ffed66'
    elif t < 5.0: return '#fff296'
    elif t < 5.5: return '#fff9d1'
    # 5.5k to 8k White
    elif t < 6.0: return '#fffced'
    elif t < 6.5: return '#f3ffed'
    elif t < 7.0: return '#ffffff'
    elif t < 7.5: return '#f4fffb'
    elif t < 8.0: return '#f4fffc'
    # 8k to 30k Light blue to dark blue
    elif t < 9.0: return '#e5fffc'
    elif t < 10: return '#c2fcf5'
    elif t < 11: return '#a0fff3'
    elif t < 12: return '#a0ebff'
    elif t < 13: return '#a0d7ff'
    elif t < 14: return '#7fc9ff'
    elif t < 15: return '#6bbdf9'
    elif t < 16: return '#49b2ff'
    elif t < 17: return '#38aaff'
    elif t < 18: return '#28a1fc'
    elif t < 19: return '#1498fc'
    elif t < 20: return '#1484fc'
    elif t < 22: return '#1465fc'
    elif t < 24: return '#144afc'
    elif t < 26: return '#142ffc'
    elif t < 28: return '#1418ff'
    elif t < 30: return '#1430ff'
    # > 30k Purple
    elif t < 32: return '#233dff'
    elif t < 35: return '#1632ff'
    elif t < 40: return '#161dff'
    elif t < 45: return '#200cff'
    elif t < 50: return '#4122f4'
    elif t < 55: return '#552ffc'
    elif t < 60: return '#6b3aff'
    else       : return '#a782ff'





#========== Class Star
#===============================================================

class Star:
  ''' Stores the magnitudes for each filter for a given star. The calibration must be done before! '''

  def SetDistance(self, d):
    ''' Sets the distance to the star, in parsec '''
    self.distance = d

  def SetMag(self, mag, filter = 'L'):
    ''' Sets the magnitude (has to be already calibrated!) for a filter '''
    if   filter == 'L': self.L = mag
    elif filter == 'R': self.R = mag
    elif filter == 'G': self.V = mag
    elif filter == 'V': self.V = mag
    elif filter == 'B': self.B = mag
    else: print '[Star.SetMag] WARNING: wrong filter (%s)'%filter

  def SetName(self, name):
    ''' Sets the name of the star... only for named stars or reference stars '''
    self.name = name

  def GetDistance(self):
    ''' Returns the distance of the star in parsec '''
    return self.distance

  def GetDistToPoint(self, x, y):
    ''' Returns the distance, in pixels, of the star to a given point in the image '''
    d2 = (x - self.GetX())*(x - self.GetX()) + (y - self.GetY())*(y - self.GetY())
    return np.sqrt(d2)

  def GetMag(self, filter = ''):
    ''' Returns the magnitude for a given filter '''
    if   filter == 'L': return self.L
    elif filter == 'R': return self.R
    elif filter == 'G': return self.V
    elif filter == 'V': return self.V
    elif filter == 'B': return self.B
    else: return self.V
    #mag = self.L if self.L != 0 else sefl.V
    #  return mag

  def GetCoor(self):
    ''' Returns the coordenates (in pixels) in the image '''
    return self.coor

  def GetName(self):
    ''' Returns the name of the star '''
    return self.name

  def GetT(self):
    ''' Returns the temperature from the B-V color index '''
    return GetTforBmV(self.GetBmV())

  def GetAbsMag(self, distance = -1):
    ''' Returns the absolute magnitude from the rel mag and the distance '''
    if distance != -1: self.SetDistance(distance)
    return GetAbsMagnitudeFromMagAndDistance(self.GetMag(), self.GetDistance())

  def GetBolometricMag(self):
    ''' Returns the bolometric magnitude from the abs magnitude and the temperature '''
    return GetBolometricMagnitude(self.GetAbsMag(), self.GetT())

  def GetLumi(self, distance = -1):
    ''' Returns the luminosity of the star from the bolometric magnitude and the distance ''' 
    if distance != -1: self.SetDistance(distance)
    return GetLumiFromMbol(self.GetBolometricMag())

  def GetRadius(self):
    ''' Returns the radius of the star using Stefan-Boltzmann law'''
    Lum = self.GetLumi()
    T = self.GetT()
    return GetStarRadius(Lum, T)

  def GetX(self):
    ''' Returns the X coordinate in the image (in pixels) '''
    return self.GetCoor()[0]

  def GetY(self):
    ''' Returns the Y coordinate in the image (in pixels) '''
    return self.GetCoor()[1]

  def GetBmV(self):
    ''' Returns the B - V color index '''
    return self.GetMag('B') - self.GetMag('V')

  def __init__(self, coor, R = -99, V = -99, B = -99, L = -99, distance = -1, name = ''):
    self.R = R
    self.V = V
    self.B = B
    self.L = L
    self.coor = coor
    self.distance = distance
    self.name = name
  
#===============================================================

class StarCluster:
  ''' Reads .FIT images, substracts bkg, gets stars in a given area, computes star objets...
      fr = StarCluster(path) # if filter can be read from the header
      fr = StarCluster(); fr.AddImage(pathToImage, 'FILTER') otherwise
      area = [x0, x1, y0, y1]
  '''

  def AddImage(self, pathToImage, filter = ''):
    ''' Add an image to the list of images of the cluster '''
    # If filter == '', takes the filter from the header
    img, f, recImg = LoadFitImage(pathToImage)
    if not filter == '': f = filter
    self.images.append(img)
    self.recImages.append(recImg)
    self.filters.append(f)
    self.removedBkg.append(False)

  def LoadImages(self):
    path = self.path
    files = os.listdir(path)
    for image in files:
      if not image.endswith('.fit') and not image.endswith('.fits') and not image.endswith('.FIT') and not image.endswith('FITS'): continue
      else:
        img, filter, recImg = LoadFitImage(path + '/' + image, self.area)
        if filter == '':
          self.info('Filter not found in file ' + image, mode = 'warn', func =  'StarCluster.LoadImages')
        elif filter == 'G': filter = 'V'
        self.images.append(img)
        self.recImages.append(recImg)
        self.filters.append(filter)
        self.removedBkg.append(False)
    self.PrintLoadedImages()
    self.info(' >> Loaded %i images!' % len(self.images), verbose = 0)

  def PrintLoadedImages(self):
    ''' Prints the filters for each loaded image '''
    for i in range(len(self.images)):
      self.info(' >> Loaded image with filter: ' + self.filters[i], verbose = 1)

  def SetArea(self, area):
    self.area = area

  def SetOutDir(self, outdir):
    self.outdir = outdir
    if not self.outdir[-1] == '/': self.outdir += '/'
    if not os.path.isdir(self.outdir): os.mkdir(self.outdir)

  def SetRefStar(self, refStar):
    self.refStar = refStar

  def AddRefStar(self, refStar):
    self.refStar.append(refStar)

  def SetClusterName(self, name):
    self.name = name

  def SetVerbose(self, v):
    self.verbose = v

  def SetDistance(self, d):
    self.distance = d

  def SetRGBtoBVR(self, val = True):
    self.doRGBtoBVR = val

  def RemoveBkgImages(self):
    i = 0
    for img in self.images:
      if not self.removedBkg[i]:
        bkg, bkg_sigma = GetBkgMedianAndSigma(img)
        self.images[i]     -= bkg
        self.removedBkg[i]  = True
      i += 1

  def SetBkgSigmaThr(self, b):
    self.bkgSigmaThr = b

  def SetSourceFWHM(self, f):
    self.fwhm = f

  def SetSourceRadius(self, s):
    self.sourceRadius = s

  def GetListOfBrightestSources(self, filter = 'V', nBrightSources = 20, area = [], options = ''):
    ''' Returns a table with the position of brigtest sources and magnitude '''
    img  = self.GetImage(filter, 'rec' in options)
    sources = GetSources(img, self.fwhm, self.bkgSigmaThr)
    X, Y, mag = GetListOfBightestSources(sources, nBrightSources, area)
    t  = '\n Brightest sources: \n'
    t += '========================== \n'
    t += '   Mag   |   X   |   Y   | \n'
    t += '-------------------------- \n'
    #    ' xxx.xxx | xxxxx | xxxxx | \n'
    for i in range(len(X)):
      t +=  ' %s | %s | %s |\n'%(SetForm(mag[i],7, '%1.2f'), SetForm(X[i], 5, '%i'), SetForm(Y[i], 5, '%i'))
    t += '========================== \n'
    self.info(t, verbose = 0)
    return t

  def DrawImage(self, filter = '', outname = '', options = '', nBrightSources = 20, area = []):
    ''' Draw the FIT with sources if 'sources' in options and the cropped image if 'rec' in options '''
    x0, x1, y0, y1 = self.GetArea()
    img  = self.GetImage(filter, 'rec' in options)

    apc = 0
    if len(self.refStar) > 0:
      xc = []; yc = []
      for s in self.refStar:
        xc.append(s.GetX())
        yc.append(s.GetY())
      posc = (xc, yc)
      apc  = CircularAperture(posc, 1.8*self.sourceRadius)
    if outname != '': fig = plt.figure(num=None, figsize=(30, 20), dpi=100, facecolor='w', edgecolor='k')
    else: fig = plt.figure()
    ax = fig.add_subplot(111)
    tr = PowerStretch(3) + PercentileInterval(98)
    img2 = tr(img)
    ax.imshow(img2, cmap='gray_r', origin='lower', interpolation='nearest', aspect='auto')
    if apc != 0:   apc.plot(color='#ff0000', lw = 2)
    if 'sources' in options:
      sources = GetSources(img, self.fwhm, self.bkgSigmaThr)
      x = sources['xcentroid']; y = sources['ycentroid']; pos = (x, y)
      bsx, bsy, mag = GetListOfBightestSources(sources, nBrightSources, area)
      bsp = (bsx, bsy)
      aper   = CircularAperture(pos, self.sourceRadius)
      aperbs = CircularAperture(bsp, 1.5*self.sourceRadius)
      if   filter == 'V':
        col = '#acef92'; col2 = '#436336'
      elif filter == 'B': 
        col = '#7ed2e5'; col2 = '#1353c1'
      elif filter == 'R': 
        col = '#ffa789'; col2 = '#c63909'
      else:
        col = '#cccc9d'; col2 = '#70705a'
      aper.plot(color = col)
      aperbs.plot(color = col2, lw = 2)
    if not 'rec' in options: 
      rect = patches.Rectangle((x0,y0),x1-x0,y1-y0,linewidth=2,edgecolor='#ffc700',facecolor='none')
      ax.add_patch(rect)
    ax.set_xlabel("X [pixels]", fontsize=10,horizontalalignment='right', position=(1,25))
    ax.set_ylabel("Y [pixels]", fontsize=10)
    if outname != '': 
      fig.savefig(self.outdir + outname+'.png')
      plt.close()
    else: plt.show()
    
  def DrawStars(self, outname = '', options = ''):
    ''' Draw the stars in positions with a color for each temperature '''
    x0, x1, y0, y1 = self.GetArea()
    x = []; y = []; T = []
    stars = self.CleanStars(self.GetSelectedStars()) if not 'rec' in options else self.CleanStars(self.GetStarsInArea(self.GetSelectedStars()))
    for s in stars:
      x.append(s.GetX())
      y.append(s.GetY())
      T.append(s.GetT())
    positions = (x, y)
    fig = plt.figure(num=None, figsize=(30, 20), dpi=100, facecolor='w', edgecolor='k')
    ax = fig.add_subplot(111)
    ax.set_xlabel("X [pixels]", fontsize=10,horizontalalignment='right', position=(1,25))
    ax.set_ylabel("Y [pixels]", fontsize=10)
    #plt.plot(x,y,'ob')
    if 'col' in options:
      ax.set_facecolor((0, 0, 0))
      ax.patch.set_facecolor((0,0,0))
      ax.spines['bottom'].set_color((1,1,1))
      ax.spines['top'].set_color((1,1,1))
      ax.spines['right'].set_color((1,1,1))
      ax.spines['left'].set_color((1,1,1))
      ax.yaxis.label.set_color((1,1,1))
      ax.xaxis.label.set_color((1,1,1))
      ax.tick_params(axis='x', colors=(1,1,1))
      ax.tick_params(axis='y', colors=(1,1,1))
      for i in range(len(x)):ax.scatter(x[i], y[i], s = 100, color = GetColor(T[i]))
      ax.scatter(x[i], y[i], s = 60, color = GetColor(T[i]))
    else:
      ax.scatter(x, y, s = 60, color = 'k')
    if not 'rec' in options:
      rect = patches.Rectangle((x0,y0),x1-x0,y1-y0,linewidth=2,edgecolor='#ffc700',facecolor='none')
      ax.add_patch(rect)
    if outname != '': 
      fig.savefig(self.outdir + outname+'.png')
      plt.close()
    else: plt.show()

  def DrawAll(self):
    filters = ['V', 'B', 'R', 'L']
    for i in range(len(filters)): 
      if not self.isThereFilter(filters[i]): filters.pop(i)
    for f in filters:
      self.DrawImage(f, outname = 'Img'+f,        options = '')
      self.DrawImage(f, outname = 'Sources'+f,    options = 'sources', nBrightSources = 20)
      self.DrawImage(f, outname = 'RecImg'+f,     options = 'rec,sources')
      self.DrawImage(f, outname = 'RecSources'+f, options = 'rec,sources', nBrightSources = 10)
    if len(self.stars) > 0:
      self.DrawStars('IdentifiedStars')
      self.DrawStars('SelectedStars', 'rec')
      self.DrawStars('IdentifiedStarsCol', 'col')
      self.DrawStars('SelectedStarsCol', 'rec,col')
    starsData    = self.PrintStarData()
    starsDataAll = self.PrintStarData(1)
    infoId = self.GetInfoIdentifiedStars()
    t = open(self.outdir + "stars_info.txt", 'w')
    t.write(infoId + '\n' + starsData)
    t.close()
    t = open(self.outdir + "stars_all_info.txt", 'w')
    t.write(infoId + '\n' + starsDataAll)
    t.close()
    self.GetCalibrationData()
    t = open(self.outdir + "calibrationData.txt", 'w')
    t.write(self.calibData)
    t.close()
    t = open(self.outdir + 'log.txt', 'w')
    t.write(self.log)
    t.close()

  def GetArea(self):
    ''' Sets the area and corrects in case of null area ''' 
    x0, x1, y0, y1 = self.area
    if not len(self.images) > 0: return self.area
    image = self.images[0]
    xmax = len(image[0])
    ymax = len(image   )
    if y1 == 0: y1 = ymax
    if x1 == 0: x1 = xmax
    return [x0, x1, y0, y1]

  def GetImage(self, filt = 'V', recImage = False):
    ''' Returns the image for the given filter ''' 
    for i in range(len(self.images)):
      img = self.images[i] if not recImage else self.recImages[i]
      fil = self.filters[i]
      if filt == fil: return img
    self.info('image not found for filter ' + filt, mode = 'warn', func = 'StarCluster.GetImage')

  def GetSources(self, filter = '', recImage = False):
    ''' Gets the sources from images '''
    img = self.GetImage(filter) if not recImage else self.GetImage(filter, 1)
    sources = GetSources(img, self.fwhm, self.bkgSigmaThr)
    return sources

  def SetStarsFilter(self, filter = 'V'):
    ''' Obtain the Star objets from one filter (V by default) '''
    self.stars = []
    if not self.isThereFilter(filter):
      self.info('filter ' + filter + ' not found!!!', mode = 'error', func = 'StarCluster.SetStarsFilters')
    sources = self.GetSources(filter)
    self.info(' >> Getting stars for filter ' + filter + '... %i stars found!!'%(len(sources)))
    for s in sources:
      x   = s['xcentroid']
      y   = s['ycentroid']
      mag = s['mag']
      a = Star(coor = [x,y])
      a.SetMag(mag, filter)
      self.stars.append(a)

  def MatchStars(self, filter):
    ''' Marches sources if a given filter with stars in a list of stars '''
    if not self.isThereFilter(filter): return
    if len(self.stars) == 0: return
    sources = self.GetSources(filter)
    self.info(' >> Getting stars for filter ' + filter + '... %i stars found!!' %len(sources))
    for s in sources:
      x   = s['xcentroid']
      y   = s['ycentroid']
      mag = s['mag']
      index = 0; indexMinD = 0; minDist = 10000
      for a in self.stars:
        dist = a.GetDistToPoint(x,y)
        if dist < minDist: 
          minDist   = dist
          indexMinD = index
        index += 1
      if minDist <= 1*self.fwhm: 
        self.stars[indexMinD].SetMag(mag, filter)

  def SetStars(self):
    ''' Get the stars from the images '''
    if len(self.stars) > 0: return
    self.SetStarsFilter('B') # Use as reference
    if self.isThereFilter('V'): self.MatchStars('V')
    if self.isThereFilter('R'): self.MatchStars('R')
    #if self.isThereFilter('L'): self.MatchStars('L')

  def GetStars(self):
    return self.stars

  def CleanStars(self, stars = ''):
    ''' Removes fake stars based on temperature and large difference of filter magnitudes '''
    cleanStars = []
    if stars == '': stars = self.stars
    lumi = []
    for s in stars:
      T = s.GetT()
      V = s.GetMag('V')
      R = s.GetMag('R')
      B = s.GetMag('B')
      if np.std([V,B,R]) > 5: continue
      if T < 2000 or T > 50000: continue
      cleanStars.append(s)
    self.info("Clean stars: %i, orig stars: %i"%(len(cleanStars),len(stars)), verbose = 0)
    return cleanStars
    
  def GetStarsInArea(self, stars = ''):
    ''' Returns the stars in the fixed area of the cluster '''
    if stars == '': stars = self.stars
    starsInArea = []
    x0, x1, y0, y1 = self.area
    for s in stars:
      x = s.GetX(); y = s.GetY()
      if x > x0 and x < x1 and y > y0 and y < y1:
        starsInArea.append(s)
    self.info('Number of stars in area: %i'%len(starsInArea), mode = 'info', func = 'StarCluster.GetStarsInArea')
    return starsInArea

  def GetSelectedStars(self, stars = ''):
    if stars == '': stars = self.stars
    selStars = []
    for s in stars:
      if s.GetMag('V') != -99 and s.GetMag('B') != -99 and s.GetMag('R') != -99:
        s.SetDistance(self.distance)
        #if self.doRGBtoBVR: ApplyRGBtoBVRtoStar(s)
        selStars.append(s)
    self.info('Number of selected stars: %i'%len(selStars), mode = 'info', func = 'StarCluster.GetSelectedStars')
    if self.doRGBtoBVR: selStars = self.CalibrateStars(selStars)
    return selStars

  def GetInfoIdentifiedStars(self):
    ''' Returns numbers of identified sources of each filter, selected stars, etc... '''
    if len(self.stars) == 0: 
      self.info('There are no stars!', func = 'StarCluster.PrintInfoStars')
      return
    nB = 0; nV = 0; nR = 0; nL = 0;
    for s in self.stars:
      if s.GetMag('V') != -99: nV += 1
      if s.GetMag('B') != -99: nB += 1
      if s.GetMag('R') != -99: nR += 1
      if s.GetMag('L') != -99: nL += 1
    t = ''
    t += '# Sources in filter V: %4i \n' %nV
    t += '# Sources in filter B: %4i \n' %nB
    t += '# Sources in filter R: %4i \n' %nR
    if nL != 0: t+= 'Sources in filter L: %4i\n' %nL
    nTotalStars  = len(self.stars)
    starsInArea = self.GetStarsInArea(self.stars)
    nStarsInArea = len(starsInArea)
    nSelStars = len(self.CleanStars(starsInArea))
    t += '\n'
    t += '# Total identified stars: %i\n' %nTotalStars
    t += '# Stars in cluster area : %i\n' %nStarsInArea
    t += '# Cleaned selected stars: %i\n' %nSelStars
    return t
    
  def GetClosestStar(self, star):
    ''' Looks for the closest star of a given star '''
    px = star.GetX()
    py = star.GetY()
    mindist = 1000000
    for s in self.stars:
      dist = s.GetDistToPoint(px, py)
      if dist < mindist:
        mindist = dist
        closestStar = s
      if mindist < self.fwhm: return s
    self.info('not found star in position [x, y] = [%i, %i]'%(px, py), func = 'StarCluster.GetStarInPos', mode = 'warn')
    return 0

  def GetCalibratingStars(self):
    ''' Returnss a list with calibrating stars '''
    cs = []
    for s in self.refStar:
      cs.append(self.GetClosestStar(s))
    return cs

  def GetMagCorrStar(self, filter, refStar = [], calibratingStar = [], kglobal = 0):
    ''' Returns constant to calibrate mangitudes from ref star for a given filter'''
    if isinstance(refStar, Star):
      mag     = refStar.GetMag(filter)
      starMag = calibratingStar.GetMag(filter)
      k = mag - starMag
    else: # It's a list of stars
      if len(refStar) == 0:         refStar = self.refStar
      if len(calibratingStar) == 0: calibratingStar = self.GetCalibratingStars()
      k = []
      for i in range(len(refStar)):
        mag     = refStar[i].GetMag(filter)
        starMag = calibratingStar[i].GetMag(filter)
        if mag == -99 or starMag == -99: continue
        k.append(mag-starMag-kglobal)
    return np.average(k), np.std(k)

  def GetGlobalMagCorr(self):
    ''' Returns a global constant '''
    filters = ['R', 'L', 'V', 'B']
    k = []
    for f in filters:
      if not self.isThereFilter(f): continue
      m,s = self.GetMagCorrStar(f)
      k.append(m)
    return np.average(k), np.std(k)

  def CalibrateStars(self, stars = 0):
    ''' Gets and applies the magnitude corrections in each filter for all the stars '''
    if stars == 0: stars = self.stars
    if self.calibData == '': self.GetCalibrationData()

    filters = ['V', 'B', 'R', 'L'];
    for i in range(len(filters)): 
      if not self.isThereFilter(filters[i]): filters.pop(i)

    calStars = self.GetCalibratingStars()

    if self.doRGBtoBVR: 
      # First, apply a global correction to all the stars
      # Then, go to BVR for the calibration stars
      # After, get the k_F calibration constant for each filter from the (BVR) calibration stars
      # Finally, for each star, apply global + go to BVR + sum k_F constant
      #ApplyRGBtoBVRtoStar(s)
      kglobal = self.GetGlobalMagCorr()[0]
      BVRcalStars = []
      for s in calStars:
        b = s.GetMag('B') + kglobal
        g = s.GetMag('V') + kglobal
        r = s.GetMag('R') + kglobal
        B = GetBGRfromRGB('B', r, g, b)
        V = GetBGRfromRGB('V', r, g, b)
        R = GetBGRfromRGB('R', r, g, b)
        BVRcalStars.append(Star([s.GetX(), s.GetY()], R=R, V=V, B=B, name=s.GetName()))
      calStars = BVRcalStars
        
    # Get calibrating constants 
    kfil = {}
    for f in filters: 
      k, sk = self.GetMagCorrStar(f, [], calStars)
      kfil[f] = k
      self.info(" >> Calibration constant = %1.2f +/- %1.2f (filter %s)"%(k,sk,f), mode='info', func='StarCluster.CalibrateStars', verbose = 2)

    for s in stars:
      if self.doRGBtoBVR: 
        kB = kfil['B']
        kV = kfil['V']
        kR = kfil['R']
        b = s.GetMag('B') + kglobal
        g = s.GetMag('V') + kglobal
        r = s.GetMag('R') + kglobal
        B = GetBGRfromRGB('B', r, g, b)
        V = GetBGRfromRGB('V', r, g, b)
        R = GetBGRfromRGB('R', r, g, b)
        s.SetMag(B+kB, 'B')
        s.SetMag(V+kV, 'V')
        s.SetMag(R+kR, 'R')
      else:
        for f in filters:
          m = s.GetMag(f)
          if m == -99: continue
          corMag = m + k
          s.SetMag(corMag, f)
    return stars
        
  def PrintStarData(self, doAllStars = False):
    ''' Returs a table with all avaliable data for the detected stars '''
    allStars = self.GetSelectedStars()
    if doAllStars: t = PrintStarData(allStars)
    else:
      stars = self.CleanStars(self.GetStarsInArea(allStars))
      t = PrintStarData(stars)
    self.info(t, verbose = 2)
    return t

  def GetCalibrationData(self):
    ''' Returs a table with the information about the mag calibration for a given ref star '''
    # Global correction
    if self.calibData != '': return self.calibData
    k = self.GetGlobalMagCorr()[0]

    # RGB corrections
    kV, sV = self.GetMagCorrStar('V', kglobal = k)
    kB, sB = self.GetMagCorrStar('B', kglobal = k)
    kR, sR = self.GetMagCorrStar('R', kglobal = k)

    # Go to BVR for calibrating stars
    calibRGB=[]
    for s in self.GetCalibratingStars():
      tG = s.GetMag('V')+k
      tR = s.GetMag('R')+k
      tB = s.GetMag('B')+k
      V = GetBGRfromRGB('V', tR, tG, tB)
      B = GetBGRfromRGB('B', tR, tG, tB)
      R = GetBGRfromRGB('R', tR, tG, tB)
      cStar = Star([s.GetX(), s.GetY()], R = R, V = V, B = B, name = s.GetName())
      calibRGB.append(cStar)
    # and re-calibrate test star
    k2V, s2V = self.GetMagCorrStar('V', [], calibRGB)
    k2B, s2B = self.GetMagCorrStar('B', [], calibRGB)
    k2R, s2R = self.GetMagCorrStar('R', [], calibRGB)

    kt = ''
    kt += '============================================\n' 
    kt += '   |       k(RGB)      |       k(BVR)      | \n'
    kt += '--------------------------------------------\n'
    k#    'B  |  xxxxx +/- xxxxx  |  xxxxx +/- xxxxx  |'
    kt += 'B  |  %s +/- %s  |  %s +/- %s  |  \n'%(SetForm(kB,5), SetForm(sB,5), SetForm(k2B,5), SetForm(s2B,5))
    kt += 'V  |  %s +/- %s  |  %s +/- %s  |  \n'%(SetForm(kV,5), SetForm(sV,5), SetForm(k2V,5), SetForm(s2V,5))
    kt += 'R  |  %s +/- %s  |  %s +/- %s  |  \n'%(SetForm(kR,5), SetForm(sR,5), SetForm(k2R,5), SetForm(s2R,5))
    kt += '============================================\n' 
    self.info(kt, verbose = 2)

    t = ''
    for refStar in self.refStar:
      name = refStar.GetName()
      calStar = self.GetClosestStar(refStar)

      # Refference
      rV = refStar.GetMag('V')
      rR = refStar.GetMag('R')
      rB = refStar.GetMag('B')

      # Raw magnitudes
      oG = calStar.GetMag('V')+k
      oR = calStar.GetMag('R')+k
      oB = calStar.GetMag('B')+k

      # Calibrated
      cG = oG+kV
      cR = oR+kR
      cB = oB+kB

      # Apply BGRfromRGB correction to test star
      o2V = GetBGRfromRGB('V', oR, oG, oB)
      o2R = GetBGRfromRGB('R', oR, oG, oB)
      o2B = GetBGRfromRGB('B', oR, oG, oB)
    
      c2V = o2V+k2V
      c2R = o2R+k2R
      c2B = o2B+k2B

      t += '\n### %s\n'%name
      t += '=======================================================================\n' 
      t += '   |  Raw mag  |  Calibrated  |  Reference  |  RGB to BVR | Calib BVR |\n'
      t += '-----------------------------------------------------------------------\n'
      #    'B  |  xxxxxxx  |  xxxxxxxxxx  |  xxxxxxxxx  |  xxxxxxxxx  |  xxxxxxx  |\n'%(oB, B, rB)
      t += 'B  |  %s  |  %s  |  %s  |  %s  |  %s  |\n'%(SetForm(oB,8), SetForm(cB,10), SetForm(rB,9), SetForm(o2B,9), SetForm(c2B,7))
      t += 'V  |  %s  |  %s  |  %s  |  %s  |  %s  |\n'%(SetForm(oG,8), SetForm(cG,10), SetForm(rV,9), SetForm(o2V,9), SetForm(c2V,7))
      t += 'R  |  %s  |  %s  |  %s  |  %s  |  %s  |\n'%(SetForm(oR,8), SetForm(cR,10), SetForm(rR,9), SetForm(o2R,9), SetForm(c2R,7))
      t += '=======================================================================\n' 

    self.calibData = kt + '\n' + t
    self.info(t, verbose = 2)
    return kt, t

  def isThereFilter(self, filter):
    for f in self.filters:
      if f == filter: return True
    return False

  # plots:
  # Phys plots, cosas corregidas y sin corregir, etc

  def __init__(self, path = '', area = [0,0,0,0], refStar = [], outdir = './', clusterName = 'cluster', distance = 0, verbose = 1, doRGBtoBVR = False):
    self.path = path
    self.area = area
    self.distance = distance
    self.images = []
    self.recImages = []
    self.filters = []
    self.removedBkg= []
    self.refStar = refStar
    self.fwhm = 5
    self.bkgSigmaThr = 3
    self.sourceRadius = 4. # pixels
    self.outdir = outdir
    if not self.outdir[-1] == '/': self.outdir += '/'
    if not os.path.isdir(self.outdir): os.mkdir(self.outdir)
    self.name = clusterName
    self.stars = []
    self.verbose = verbose
    self.log = ''
    self.doRGBtoBVR = doRGBtoBVR
    self.calibData = ''
    if self.path != '': self.LoadImages()

  def info(self, t, verbose = 1, mode = '', func = ''):
    pre = ''
    if func != '': pre = '[' + func + '] '
    if   mode == 'warn' : pre += 'WARNING: '
    elif mode == 'error': pre += 'ERROR: '
    elif mode == 'info' : pre += 'INFO: '
    text = pre + t
    if   mode == 'warn' or mode == 'error': print text
    elif self.verbose >= verbose:           print text
    self.log += text + '\n'

