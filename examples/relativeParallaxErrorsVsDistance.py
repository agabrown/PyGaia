#!/usr/bin/env python

# Plot relative parallax errors as a function of distance for stars of a given spectral type.

import numpy as np
import matplotlib.pyplot as plt
from collections import OrderedDict
from pygaia.photometry.transformations import gminvFromVmini
from pygaia.errors.astrometric import parallaxErrorSkyAvg

from os import environ as env
from matplotlib import rc
try:
  import argparse
except ImportError:
  raise ImportError("""The argparse module does not seem to be available.\n Either upgrade to python 2.7
  or adapt the script to using the optparse module instead.""")

# Configure matplotlib
rc('text', usetex=True)
rc('font', family='serif', size=16)
rc('xtick.major', size='12')
rc('xtick.minor', size='6')
rc('ytick.major', size='12')
rc('ytick.minor', size='6')
rc('lines', linewidth=2)
rc('axes', linewidth=2)

def makePlot(pdf=False, png=False):
  """
  Plot relative parallax errors as a function of distance for stars of a given spectral type.

  Parameters
  ----------

  args - command line arguments
  """
  logdistancekpc = np.linspace(-1,np.log10(20.0),100)
  sptVabsAndVmini=OrderedDict([('K0V',(5.58,0.87)), ('G5V',(4.78,0.74)), ('G0V',(4.24,0.67)),
    ('F5V',(3.50,0.50)), ('F0V',(2.98,0.38)), ('RC',(0.8,1.0))])
  lines={}

  fig=plt.figure(figsize=(10,6.5))
  currentAxis=plt.gca()

  for spt in sptVabsAndVmini.keys():
    vmag=sptVabsAndVmini[spt][0]+5.0*logdistancekpc+10.0
    indices=(vmag>14) & (vmag<16)
    gmag=vmag+gminvFromVmini(sptVabsAndVmini[spt][1])
    parerrors=parallaxErrorSkyAvg(gmag,sptVabsAndVmini[spt][1])
    relparerrors=parerrors*10**logdistancekpc/1000.0
    plt.loglog(10**logdistancekpc, relparerrors,'--k',lw=1)
    plt.loglog(10**logdistancekpc[indices], relparerrors[indices],'-',label=spt)
  plt.xlim(0.1,20.0)
  plt.ylim(0.001,0.5)
  plt.text(0.9, 0.05,'Colours indicate $14<V<16$',
     horizontalalignment='right',
     verticalalignment='bottom',
     transform = currentAxis.transAxes)
  plt.legend(loc=2)
  plt.xlabel('distance [kpc]')
  plt.ylabel('$\\sigma_\\varpi/\\varpi$')
  plt.grid(which='both')
  
  if (args['pdfOutput']):
    plt.savefig('RelativeParallaxErrorsVsDist.pdf')
  elif (args['pngOutput']):
    plt.savefig('RelativeParallaxErrorsVsDist.png')
  else:
    plt.show()

def parseCommandLineArguments():
  """
  Set up command line parsing.
  """
  parser = argparse.ArgumentParser(description="""Plot lines of constant distance for stars a given G
  magnitude in the Mv vs (V-I) diagram""")
  parser.add_argument("-p", action="store_true", dest="pdfOutput", help="Make PDF plot")
  parser.add_argument("-b", action="store_true", dest="pngOutput", help="Make PNG plot")
  args=vars(parser.parse_args())
  return args

if __name__ in ('__main__'):
  args=parseCommandLineArguments()
  makePlot(pdf=args['pdfOutput'], png=args['pngOutput'])
