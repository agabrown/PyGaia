#!/usr/bin/env python

# Plot photometric errors as a function of V according to the interpolation formulae on the Gaia
# science performance web-pages.

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from pygaia.errors.photometric import gMagnitudeError, bpMagnitudeError, rpMagnitudeError
from pygaia.errors.photometric import gMagnitudeErrorEoM, bpMagnitudeErrorEoM, rpMagnitudeErrorEoM
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
  Make the plot with photometry performance predictions.

  :argument args: command line arguments
  """
  gmag=np.linspace(3.0,20.0,171)

  vmini = args['vmini']
  
  vmag=gmag-gminvFromVmini(vmini)
  
  if args['eom']:
      sigmaG = gMagnitudeErrorEoM(gmag)
      sigmaGBp = bpMagnitudeErrorEoM(gmag, vmini)
      sigmaGRp = rpMagnitudeErrorEoM(gmag, vmini)
      yminmax = (1.0-4,0.1)
  else:
      sigmaG = gMagnitudeError(gmag)
      sigmaGBp = bpMagnitudeError(gmag, vmini)
      sigmaGRp = rpMagnitudeError(gmag, vmini)
      yminmax = (1.0-4,1)

  fig=plt.figure(figsize=(10,6.5))
  
  if (args['vmagAbscissa']):
    plt.semilogy(vmag, sigmaG, 'k', label='$\\sigma_G$')
    plt.semilogy(vmag, sigmaGBp, 'b', label='$\\sigma_{G_\\mathrm{BP}}$'+' for $(V-I)={0}$'.format(vmini))
    plt.semilogy(vmag, sigmaGRp, 'r', label='$\\sigma_{G_\\mathrm{RP}}$'+' for $(V-I)={0}$'.format(vmini))
    plt.xlim((6,20))
    #plt.ylim(yminmax)
    plt.legend(loc=0)
    plt.xlabel('$V$ [mag]')
  else:
    ax=fig.add_subplot(111)
    plt.semilogy(gmag, sigmaG, 'k', label='$\\sigma_G$')
    plt.semilogy(gmag, sigmaGBp, 'b', label='$\\sigma_{G_\\mathrm{BP}}$'+' for $(V-I)={0}$'.format(vmini))
    plt.semilogy(gmag, sigmaGRp, 'r', label='$\\sigma_{G_\\mathrm{RP}}$'+' for $(V-I)={0}$'.format(vmini))
    plt.xlim((6,20))
    #plt.ylim(yminmax)
    plt.legend(loc=0)
    plt.xlabel('$G$ [mag]')
  
  plt.xticks(np.arange(6,20,2))
  ax = plt.gca().yaxis 
  #ax.set_major_formatter(matplotlib.ticker.ScalarFormatter())
  #plt.ticklabel_format(axis='y',style='plain')
  plt.grid(which='both')
  plt.ylabel('Photometric error [mag]')
  if args['eom']:
      plt.title('End-of-mission mean photometry: sky averaged errors for $(V-I)={0}$'.format(vmini), fontsize=14)
  else:
      plt.title('Single-FoV-transit photometry: sky averaged errors for $(V-I)={0}$'.format(vmini), fontsize=14)
  
  basename = 'PhotometricErrors'
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
  parser = argparse.ArgumentParser(description="Plot predicted single-transit Gaia sky averaged photometric errors as a function of V or G")
  parser.add_argument("vmini", help="""(V-I) colour of source""", type=float)
  parser.add_argument("-e", action="store_true", dest="eom", help="Plot end-of-mission errors")
  parser.add_argument("-p", action="store_true", dest="pdfOutput", help="Make PDF plot")
  parser.add_argument("-b", action="store_true", dest="pngOutput", help="Make PNG plot")
  parser.add_argument("-v", action="store_true", dest="vmagAbscissa", help="Plot performance vs V instead of G")
  args=vars(parser.parse_args())
  return args

if __name__ in ('__main__'):
  args=parseCommandLineArguments()
  makePlot(args)
