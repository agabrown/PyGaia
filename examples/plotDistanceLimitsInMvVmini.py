#!/usr/bin/env python

# Plot the lines of constant distance for stars at G=20 in a colour magnitude diagram (Mv vs V-I).

import numpy as np
import matplotlib.pyplot as plt
from pygaia.photometry.transformations import gminvFromVmini, vminGrvsFromVmini

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

def makePlot(gmag, pdf=False, png=False, rvs=False):
  """
  Make a plot of a Mv vs (V-I) colour magnitude diagram containing lines of constant distance for stars
  at G=20. This will give an idea of the reach of Gaia.

  Parameters
  ----------

  args - command line arguments
  """
  vmini = np.linspace(-0.5,4.0,100)
  if (rvs):
    gminv = -vminGrvsFromVmini(vmini)
  else:
    gminv = gminvFromVmini(vmini)
  mvlimit100pc = gmag-5.0*np.log10(100.0)+5.0-gminv
  mvlimit1kpc = gmag-5.0*np.log10(1000.0)+5.0-gminv
  mvlimit10kpc = gmag-5.0*np.log10(10000.0)+5.0-gminv

  fig=plt.figure(figsize=(8,8))
  plt.plot(vmini,mvlimit100pc,'b')
  plt.text(vmini[50]-0.4,mvlimit100pc[50],"$d=100$ pc", horizontalalignment='right', va='top')
  plt.plot(vmini,mvlimit1kpc,'r')
  plt.text(vmini[50]-0.4,mvlimit1kpc[50],"$d=1000$ pc", horizontalalignment='right', va='top')
  plt.plot(vmini,mvlimit10kpc,'g')
  plt.text(vmini[50]-0.4,mvlimit10kpc[50],"$d=10000$ pc", horizontalalignment='right', va='top')
  ax=plt.gca()
  ax.set_ylim(ax.get_ylim()[::-1])
  plt.xlabel("$(V-I)$")
  plt.ylabel("$M_V$")
  if (rvs):
    plt.title("Distance limits for $G_\\mathrm{RVS}"+"={0}$".format(gmag))
  else:
    plt.title("Distance limits for $G={0}$".format(gmag))
  
  if (args['pdfOutput']):
    plt.savefig('GaiaSurveyLimits.pdf')
  elif (args['pngOutput']):
    plt.savefig('GaiaSurveyLimits.png')
  else:
    plt.show()

def parseCommandLineArguments():
  """
  Set up command line parsing.
  """
  parser = argparse.ArgumentParser(description="""Plot lines of constant distance for stars a given G
  magnitude in the Mv vs (V-I) diagram""")
  parser.add_argument("gmag", help="G-band magnitude", type=float)
  parser.add_argument("-p", action="store_true", dest="pdfOutput", help="Make PDF plot")
  parser.add_argument("-b", action="store_true", dest="pngOutput", help="Make PNG plot")
  parser.add_argument("-r", action="store_true", dest="useGrvs", help="Magnitude limit is in Grvs")
  args=vars(parser.parse_args())
  return args

if __name__ in ('__main__'):
  args=parseCommandLineArguments()
  makePlot(args['gmag'], pdf=args['pdfOutput'], png=args['pngOutput'], rvs=args['useGrvs'])
