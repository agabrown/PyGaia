#!/usr/bin/env python

# Plot radial velocity errors as a function of V according to the interpolation
# formulae on the Gaia science performance web-pages.

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from pygaia.errors.spectroscopic import vradErrorSkyAvg
from pygaia.photometry.transformations import vminGrvsFromVmini
from pygaia.photometry.utils import vminiFromSpt

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
rc('font', family='serif', size=16)
rc('xtick.major', size='12')
rc('xtick.minor', size='6')
rc('ytick.major', size='12')
rc('ytick.minor', size='6')
rc('lines', linewidth=2)
rc('axes', linewidth=2)

def makePlot(args):
  """
  Make the plot with radial velocity performance predictions.

  :argument args: command line arguments
  """
  gRvs=np.linspace(5.7,16.1,101)

  spts=['B0V', 'B5V', 'A0V', 'A5V', 'F0V', 'G0V',
        'G5V', 'K0V', 'K1IIIMP', 'K4V', 'K1III']

  fig=plt.figure(figsize=(10,6.5))
  deltaHue = 240.0/(len(spts)-1)
  hsv=np.zeros((1,1,3))
  hsv[0,0,1]=1.0
  hsv[0,0,2]=0.9
  count=0
  for spt in spts:
    hsv[0,0,0]=(240-count*deltaHue)/360.0
    vmag = vminGrvsFromVmini(vminiFromSpt(spt)) + gRvs
    vradErrors = vradErrorSkyAvg(vmag, spt)
    plt.plot(vmag, vradErrors, '-', label=spt, color=hsv_to_rgb(hsv)[0,0,:])
    count+=1
  plt.grid(which='both')
  plt.xlim(9,17.5)
  plt.ylim(0,20)
  plt.xticks(np.arange(9,18,1))
  plt.yticks(np.arange(0,20.5,5))
  plt.xlabel('$V$ [mag]')
  plt.ylabel('End-of-mission radial velocity error [km s$^{-1}$]')
  leg=plt.legend(loc=0,  handlelength=2.0, labelspacing=0.10)
  for t in leg.get_texts():
    t.set_fontsize(12)

  if (args['pdfOutput']):
    plt.savefig('RadialVelocityErrors.pdf')
  elif (args['pngOutput']):
    plt.savefig('RadialVelocityErrors.png')
  else:
    plt.show()

def parseCommandLineArguments():
  """
  Set up command line parsing.
  """
  parser = argparse.ArgumentParser(description="Plot predicted Gaia sky averaged radial velocity errors as a function of V")
  parser.add_argument("-p", action="store_true", dest="pdfOutput", help="Make PDF plot")
  parser.add_argument("-b", action="store_true", dest="pngOutput", help="Make PNG plot")
  parser.add_argument("-g", action="store_true", dest="grvsAbscissa", help="Plot performance vs Grvs instead of V")
  args=vars(parser.parse_args())
  return args

if __name__ in ('__main__'):
  args=parseCommandLineArguments()
  makePlot(args)
