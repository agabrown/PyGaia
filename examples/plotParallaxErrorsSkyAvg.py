#!/usr/bin/env python

# Plot parallax errors as a function of V according to the interpolation
# formulae on the Gaia science performance web-pages.

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from pygaia.errors.astrometric import parallaxErrorSkyAvg, parallaxErrorSkyAvgAltStartGate
from pygaia.photometry.utils import vminiFromSpt
from pygaia.photometry.transformations import gminvFromVmini

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

def makePlot(args):
  """
  Make the plot with parallax performance predictions.

  :argument args: command line arguments
  """
  gmag=np.linspace(5.7,20.0,101)

  vminiB1V=vminiFromSpt('B1V')
  vminiG2V=vminiFromSpt('G2V')
  vminiM6V=vminiFromSpt('M6V')
  
  vmagB1V=gmag-gminvFromVmini(vminiB1V)
  vmagG2V=gmag-gminvFromVmini(vminiG2V)
  vmagM6V=gmag-gminvFromVmini(vminiM6V)
  
  sigparB1V=parallaxErrorSkyAvg(gmag,vminiB1V)
  sigparB1VBrightBox=1.156*sigparB1V
  indices=(sigparB1VBrightBox<14.0)
  sigparB1VBrightBox[indices]=14.0
  sigparB1VAltStartGate=parallaxErrorSkyAvgAltStartGate(gmag,vminiB1V)
  
  sigparG2V=parallaxErrorSkyAvg(gmag,vminiG2V)
  sigparG2VBrightBox=1.156*sigparG2V
  indices=(sigparG2VBrightBox<14.0)
  sigparG2VBrightBox[indices]=14.0
  sigparG2VAltStartGate=parallaxErrorSkyAvgAltStartGate(gmag,vminiG2V)
  
  sigparM6V=parallaxErrorSkyAvg(gmag,vminiM6V)
  sigparM6VBrightBox=1.156*sigparM6V
  indices=(sigparM6VBrightBox<14.0)
  sigparM6VBrightBox[indices]=14.0
  sigparM6VAltStartGate=parallaxErrorSkyAvgAltStartGate(gmag,vminiM6V)
  
  fig=plt.figure(figsize=(10,6.5))
  
  if (args['gmagAbscissa']):
    plt.semilogy(gmag, sigparB1V, 'b', label='B1V')
    plt.semilogy(gmag, sigparG2V, 'g', label='G2V')
    plt.semilogy(gmag, sigparM6V, 'r', label='M6V')
    plt.xlim((5,20))
    plt.ylim((4,400))
  else:
    ax=fig.add_subplot(111)
    plt.semilogy(vmagB1V, sigparB1V, 'b', label='B1V')
    #plt.semilogy(vmagG2V, sigparG2V, 'g', label='G2V')
    plt.semilogy(vmagM6V, sigparM6V, 'r', label='M6V')
    plt.fill_between(vmagB1V, 0.7*sigparB1V, sigparB1VBrightBox, color='b', alpha=0.3)
    plt.fill_between(vmagM6V, 0.7*sigparM6V, sigparM6VBrightBox, color='r', alpha=0.3)
    plt.xlim((5,22.5))
    plt.ylim((4,1000))
    plt.text(17.5,190,'B1V',color='b')
    plt.text(18,20,'M6V',color='r')
    plt.text(7,17,'calibration noise floor', size=12, bbox=dict(boxstyle="round,pad=0.3",
                       ec=(0.0, 0.0, 0.0),
                       fc=(1.0, 1.0, 1.0),
                       ))
    plt.text(14.75,80,'photon noise', rotation=45, size=12, bbox=dict(boxstyle="round,pad=0.3",
                       ec=(0.0, 0.0, 0.0),
                       fc=(1.0, 1.0, 1.0),
                       ))
    ax.annotate('non-uniformity\nover the sky', xy=(21.5, 160),  xycoords='data',
                  xytext=(21.5,50), textcoords='data', ha='center', size='12',
                  bbox=dict(boxstyle="round,pad=0.3",ec=(0,0,0),fc=(1,1,1)),
                  arrowprops=dict(facecolor='black', shrink=0.15, width=1,
                    headwidth=6),
                  horizontalalignment='right', verticalalignment='top',
                  )
    ax.annotate('', xy=(21.5, 250),  xycoords='data',
                  xytext=(21.5,800), textcoords='data', ha='center', size='12',
                  arrowprops=dict(facecolor='black', shrink=0.15, width=1,
                    headwidth=6),
                  horizontalalignment='right', verticalalignment='bottom',
                  )
  
  plt.xticks(np.arange(6,24,2))
  ax = plt.gca().yaxis 
  ax.set_major_formatter(matplotlib.ticker.ScalarFormatter())
  plt.ticklabel_format(axis='y',style='plain')
  plt.grid(which='both')
  plt.xlabel('$V$ [mag]')
  plt.ylabel('End-of-mission parallax standard error [$\mu$as]')
  
  if (args['pdfOutput']):
    plt.savefig('ParallaxErrors.pdf')
  elif (args['pngOutput']):
    plt.savefig('ParallaxErrors.png')
  else:
    plt.show()

def parseCommandLineArguments():
  """
  Set up command line parsing.
  """
  parser = argparse.ArgumentParser(description="Plot predicted Gaia sky averaged parallax errors as a function of V")
  parser.add_argument("-p", action="store_true", dest="pdfOutput", help="Make PDF plot")
  parser.add_argument("-b", action="store_true", dest="pngOutput", help="Make PNG plot")
  parser.add_argument("-g", action="store_true", dest="gmagAbscissa", help="Plot performance vs G instead of V")
  args=vars(parser.parse_args())
  return args

if __name__ in ('__main__'):
  args=parseCommandLineArguments()
  makePlot(args)
