'''
  @ autor: Xuan
  @ mail: gonzalezrodrigo@uniovi.es

  Lee datos de estrellas de una tabla y produce bonitos diagramas
  color-magnitud, Lumi-temperatura, etc.

'''
from StarCluster import GetColor, GetTforBmV
import numpy as np
import os, sys
import matplotlib.pyplot as plt
import matplotlib.patches as patches

class StarPlotter:

  ################################# 'Set' methods

  def SetOutDir(self, outdir):
    self.outdir = outdir
    if not self.outdir[-1] == '/': outdir += '/'

  def SetFigSize(self, x, y):
    self.figsize = (x,y)

  def SetSecondYaxis(self, val = True):
    self.doTaxis = val

  def SetNticksX(self, nX):
    self.nTicksX = nX

  def SetNticksY(self, nY):
    self.nTicksY = nY

  def SetBaseStarArea(self, b):
    self.BaseStarArea = b

  ################################# 'Get' methods

  def GetStarIndices(self):
    return self.index

  def GetNstars(self):
    return len(self.indes)

  def GetX(self):
    return self.X

  def GetY(self):
    return self.Y

  def GetCoor(self):
    return [self.X, self.Y]

  def GetMag(self, f = ''):
    if f == '': f = 'V'
    if   f == 'V': return self.V
    elif f == 'R': return self.R
    elif f == 'B': return self.B
    else:
      print  'WARNING: filter \'%s\' not found...'%f
      return []

  def GetAbsMag(self):
    return self.absMag

  def GetBolMag(self):
    return self.BolMag

  def GetLumi(self):
    return self.Lumi

  def GetBmV(self):
    return self.BmV

  def GetT(self):
    return self.T

  def GetStarRadius(self):
    return self.Rad

  #################################################

  def GetData(self):
    ''' Autimatically reads the data; a table with all the data must exists... run StarCluster first!! '''
    f = open(self.path)
    t = f.readlines()
    val = []
    for line in t:
      if line.startswith('='): continue
      if line.startswith('-'): continue
      if line.startswith('#'): continue
      if line.replace(' ', '').startswith('id'): continue
      l = line.replace(' ', '').split('|')[:-1]
      if len(l) != 0: val.append(l)
    val = np.transpose(val)
    self.index  = [float(i) for i in val[0]]
    self.X      = [float(i) for i in val[1]]
    self.Y      = [float(i) for i in val[2]]
    self.B      = [float(i) for i in val[3]]
    self.V      = [float(i) for i in val[4]]
    self.R      = [float(i) for i in val[5]]
    self.absMag = [float(i) for i in val[6]]
    self.BolMag = [float(i) for i in val[7]]
    self.Lumi   = [float(i) for i in val[8]]
    self.BmV    = [float(i) for i in val[9]]
    self.T      = [float(i) for i in val[10]]
    self.Rad    = [float(i) for i in val[11]]


  def SetPath(self, path):
    ''' Set the path to the dir where the .txt file is '''
    if os.path.isdir(path):
      fCand  = []; candLines = []; goodFile = ''
      if not path[-1] == '/': path += '/'
      # Selecting candidates
      for f in os.listdir(path):
        if f.endswith('.txt'): fCand.append(f)
      # Check if contains star data... if so: count number of lines
      for table in fCand:
        nLines = 99999999999999
        with open(path + table, 'r') as f:
          lines = f.readlines()[0:100]
          for l in lines:
            if 'id' in l and 'Lumi' in l and 'B-V' in l: 
              nLines = len(lines)
              break
        candLines.append(nLines)
      nMinLines = candLines[0]; goodFile = fCand[0]
      for i in range(len(candLines)):
        if candLines[i] < nMinLines:
          nMinLines = candLines[i]      
          goodFile = fCand[i]
      if goodFile == '':
        print 'ERROR: not found valid file in path: %s'%path
      else:
        print ' >> Data from file: %s'%(path+goodFile)
        return path+goodFile
    else:
      if os.path.isfile(path): return path
      if os.path.isfile(path + '.txt'): return path + '.txt'

  def CreateHRdiagram(self, outname = '', options = ''):
    ''' Creates a nice HR diagram '''
    self.CreateColorMagDiagram(outname, 'lum,'+options)

  def CreateColorMagDiagram(self, outname = '', options = ''):
    ''' Creates a nice magnitude-color diagram '''
    lumi = self.Lumi; colo = self.BmV; vmag = self.V; T = self.T
    x = colo;
    y = vmag if not 'lum' in options else lumi
    fig = plt.figure(num=None, figsize=self.figsize, dpi=100, facecolor='w', edgecolor='k')

    # Create another X axis for temperatures
    ax1 = fig.add_subplot(111)
    #ax1.set_xticks(np.array([0,0.5,1]))

    # Invert the Y axis (the magnitude, as it's inverse to the luminosity)
    if not 'lum' in options: plt.gca().invert_yaxis()
    #ax1.set_xlim([-0.34,0.9247])

    # Set labels to axes
    ax1.set_xlabel("B-V", fontsize=20,horizontalalignment='right', position=(1,25))
    ax1.set_ylabel("Magnitude",fontsize=20)
    if 'lum' in options:
      ax1.set_ylabel("Luminosity",fontsize=20)
      ax1.set_yscale('log')

    # Set colors
    ax1.set_facecolor((0, 0, 0))
    fig.patch.set_facecolor((0,0,0))
    ax1.spines['bottom'].set_color((1,1,1))
    ax1.spines['top'].set_color((1,1,1))
    ax1.spines['right'].set_color((1,1,1))
    ax1.spines['left'].set_color((1,1,1))
    ax1.yaxis.label.set_color((1,1,1))
    ax1.xaxis.label.set_color((1,1,1))
    ax1.tick_params(axis='x', colors=(1,1,1))
    ax1.tick_params(axis='y', colors=(1,1,1))

    if self.doTaxis:
      ax2 = ax1.twiny()
      #Tticks = np.array([-0.339199, -0.0294, 0.3434, 0.9247, 1.791])
      # Set the second X axis
      #ax2.set_xlim(ax1.get_xlim())
      #ax2.set_xticks(Tticks)
      #mn, mx = ax2.get_xlim()
      #ax2.set_xlim(GetT(mn), GetT(mx))
      ax2.set_xlabel("Temperature (K)", fontsize = 20, horizontalalignment='right', position=(1,25))
      ax2.spines['bottom'].set_color((1,1,1))
      ax2.spines['top'].set_color((1,1,1))
      ax2.spines['right'].set_color((1,1,1))
      ax2.spines['left'].set_color((1,1,1))
      ax2.xaxis.label.set_color((1,1,1))
      ax2.tick_params(axis='x', colors=(1,1,1))
      #ax2.get_yaxis().set_label_coords()

    for i in range(len(x)):
      ax1.scatter(x[i], y[i], s = self.BaseStarArea, color = GetColor(T[i]))

    ymin, ymax = ax1.get_ylim()
    xmin, xmax = ax1.get_xlim()
    if not 'lum' in options: ax1.set_yticks(np.round(np.linspace(ymin, ymax, self.nTicksY), 2))
    ax1.set_xticks(np.round(np.linspace(xmin, xmax, self.nTicksX), 2))
    if self.doTaxis:
      ax2.set_xticks([GetTforBmV(i) for i in np.round(np.linspace(xmin, xmax, self.nTicksX/2), 2)])
      ax2.invert_xaxis()

    # Show or save
    if outname != '': 
      fig.savefig(self.outdir + outname+'.png', facecolor=fig.get_facecolor(), edgecolor='none')
      plt.close()
    else: plt.show()

  def DrawRmVtoBmV(self, outname = ''):
    ''' Draw the R-V to B-V plot '''
    fig = plt.figure(num=None, figsize=self.figsize, dpi=100, facecolor='w', edgecolor='k')
    x = [a-b for (a,b) in zip(self.B,self.V)]
    y = [a-b for (a,b) in zip(self.V,self.R)]
    plt.plot(x,y, 'ob')

    if outname != '': 
      fig.savefig(self.outdir + outname+'.png', facecolor=fig.get_facecolor(), edgecolor='none')
      plt.close()
    else: plt.show()
    
  def __init__(self, pathToTable, outdir = './', addTaxis = True, figsize = (30, 20), baseStarArea = 100):
    self.path = self.SetPath(pathToTable)
    self.GetData()
    self.figsize = figsize
    self.doTaxis = addTaxis
    self.nTicksX = 12
    self.nTicksY = 12
    self.BaseStarArea = baseStarArea
    self.outdir = outdir
    if not self.outdir[-1] == '/': outdir += '/'

