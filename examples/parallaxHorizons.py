"""
Plot the parallax horizons for stars of various spectral types. The parallax horizon defines out to what
distance the star can be seen for a given relative parallax accuracy.

Anthony Brown 2013
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from pygaia.errors.astrometric import parallaxErrorSkyAvg
from pygaia.photometry.utils import vminiFromSpt, gabsFromSpt, vabsFromSpt
from pygaia.photometry.transformations import gminvFromVmini

from os import environ as env
from matplotlib import rc
from matplotlib.colors import hsv_to_rgb

try:
  import argparse
except ImportError:
  raise ImportError("""The argparse module does not seem to be available.\n Either upgrade to python 2.7
  or adapt the script to using the optparse module instead.""")

# Configure matplotlib
rc('text', usetex=True)
rc('font', family='serif', size=18)
rc('xtick.major', size='12')
rc('xtick.minor', size='6')
rc('ytick.major', size='12')
rc('ytick.minor', size='6')
rc('lines', linewidth=2)
rc('axes', linewidth=2)

def makePlot(args):
  """
  Make the plot with parallax horizons. The plot shows V-band magnitude vs distance for a number of
  spectral types and over the range 5.7<G<20. In addition a set of crudely drawn contours show the points
  where 0.1, 1, and 10 per cent relative parallax accracy are reached.

  Parameters
  ----------
  
  args - Command line arguments.
  """
  distances = 10.0**np.linspace(1,6,10001)
  av = args['extinction']
  ai = 0.479*av #Cardelli et al R=3.1

  spts = ['B0I', 'B1V', 'G2V', 'K4V', 'M0V', 'M6V', 'K1III', 'M0III']
  pointOnePercD = []
  pointOnePercV = []
  onePercD = []
  onePercV = []
  tenPercD = []
  tenPercV = []
  vabsPointOnePerc = []
  vabsOnePerc = []
  vabsTenPerc = []

  fig=plt.figure(figsize=(11,7.8))
  deltaHue = 240.0/(len(spts)-1)
  hues = (240.0-np.arange(len(spts))*deltaHue)/360.0
  hsv=np.zeros((1,1,3))
  hsv[0,0,1]=1.0
  hsv[0,0,2]=0.9
  for hue,spt in zip(hues, spts):
    hsv[0,0,0]=hue
    vmags = vabsFromSpt(spt)+5.0*np.log10(distances)-5.0+av
    vmini=vminiFromSpt(spt)+av-ai
    #gmags = gabsFromSpt(spt)+5.0*np.log10(distances)-5.0
    gmags = vmags + gminvFromVmini(vmini)
    relParErr = parallaxErrorSkyAvg(gmags,vmini)*distances/1.0e6
    observed = (gmags>=5.7) & (gmags<=20.0)
    relParErrObs = relParErr[observed]
    # Identify the points where the relative parallax accuracy is 0.1, 1, or 10 per cent.
    if (relParErrObs.min()<0.001):
      index = len(relParErrObs[relParErrObs<=0.001])-1
      pointOnePercD.append(distances[observed][index])
      pointOnePercV.append(vmags[observed][index])
      vabsPointOnePerc.append(vabsFromSpt(spt))
    if (relParErrObs.min()<0.01):
      index = len(relParErrObs[relParErrObs<=0.01])-1
      onePercD.append(distances[observed][index])
      onePercV.append(vmags[observed][index])
      vabsOnePerc.append(vabsFromSpt(spt))
    if (relParErrObs.min()<0.1):
      index = len(relParErrObs[relParErrObs<=0.1])-1
      tenPercD.append(distances[observed][index])
      tenPercV.append(vmags[observed][index])
      vabsTenPerc.append(vabsFromSpt(spt))
    plt.semilogx(distances[observed], vmags[observed], '-', label=spt, color=hsv_to_rgb(hsv)[0,0,:])
    if (spt=='B0I'):
      plt.text(distances[observed][-1]-1.0e5, vmags[observed][-1], spt, horizontalalignment='right',
          verticalalignment='bottom', fontsize=14)
    else:
      plt.text(distances[observed][-1], vmags[observed][-1], spt, horizontalalignment='center',
          verticalalignment='bottom', fontsize=14)

  # Draw the "contours" of constant relative parallax accuracy.
  pointOnePercD = np.array(pointOnePercD)
  pointOnePercV = np.array(pointOnePercV)
  indices = np.argsort(vabsPointOnePerc)
  plt.semilogx(pointOnePercD[indices],pointOnePercV[indices],'k--')
  plt.text(pointOnePercD[indices][-1]*1.2,pointOnePercV[indices][-1]-2.5,"$0.1$\\%", ha='right', size=16,
      bbox=dict(boxstyle="round, pad=0.3", ec=(0.0, 0.0, 0.0), fc=(1.0, 1.0, 1.0),))

  onePercD = np.array(onePercD)
  onePercV = np.array(onePercV)
  indices = np.argsort(vabsOnePerc)
  plt.semilogx(onePercD[indices],onePercV[indices],'k--')
  plt.text(onePercD[indices][-1]*1.2,onePercV[indices][-1]-2.5,"$1$\\%", ha='right', size=16,
      bbox=dict(boxstyle="round, pad=0.3", ec=(0.0, 0.0, 0.0), fc=(1.0, 1.0, 1.0),))

  tenPercD = np.array(tenPercD)
  tenPercV = np.array(tenPercV)
  indices = np.argsort(vabsTenPerc)
  plt.semilogx(tenPercD[indices],tenPercV[indices],'k--')
  plt.text(tenPercD[indices][-1]*1.5,tenPercV[indices][-1]-2.5,"$10$\\%", ha='right', size=16,
      bbox=dict(boxstyle="round, pad=0.3", ec=(0.0, 0.0, 0.0), fc=(1.0, 1.0, 1.0),))

  plt.title('Parallax relative accuracy horizons ($A_V={0}$)'.format(av))

  plt.xlabel('Distance [pc]')
  plt.ylabel('V')
  plt.grid()
  #leg=plt.legend(loc=4, fontsize=14, labelspacing=0.5)
  plt.ylim(5,26)
  
  basename='ParallaxHorizons'
  if (args['pdfOutput']):
    plt.savefig(basename+'.pdf')
  elif (args['pngOutput']):
    plt.savefig(basename+'.png')
  else:
    plt.show()

def parseCommandLineArguments():
  """
  Set up command line parsing.
  """
  parser = argparse.ArgumentParser(description="Plot parallax horizons for various spectral types.")
  parser.add_argument("-a", dest="extinction", type=float, help="Value of extinction Av", default=0.0)
  parser.add_argument("-p", action="store_true", dest="pdfOutput", help="Make PDF plot")
  parser.add_argument("-b", action="store_true", dest="pngOutput", help="Make PNG plot")
  args=vars(parser.parse_args())
  return args

if __name__ in ('__main__'):
  args=parseCommandLineArguments()
  makePlot(args)
