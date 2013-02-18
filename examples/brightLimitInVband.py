"""
For a given G-band bright limit for Gaia detection/confirmation, plot the bright limit in the V-band as a
function of (V-I).
"""

from pygaia.photometry.transformations import gminvFromVmini
import numpy as np
import matplotlib.pyplot as plt
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

def plotBrightLimitInV(gBright, pdf=False, png=False):
  """
  Plot the bright limit of Gaia in V as a function of (V-I).

  Parameters
  ----------

  gBright - The bright limit of Gaia in G
  """
  vmini=np.linspace(0.0,6.0,1001)
  gminv=gminvFromVmini(vmini)
  vBright=gBright-gminv

  fig=plt.figure(figsize=(10,6.5))
  plt.plot(vmini,vBright,'b-')
  plt.xlabel('$(V-I)$')
  plt.ylabel('Bright limit of Gaia in $V$')
  plt.xlim(0,6)
  plt.ylim(5,11)
  plt.grid(which='both')
  plt.title("Bright limit in $G$: {0}".format(gBright))

  if (pdf):
    plt.savefig('VBandBrightLimit.pdf')
  elif (png):
    plt.savefig('VBandBrightLimit.png')
  else:
    plt.show()

def parseCommandLineArguments():
  """
  Set up command line parsing.
  """
  parser = argparse.ArgumentParser(description="""Plot the bright limit in the V-band, as a function of
  (V-I) for a given bright limit in G""")
  parser.add_argument("gmag", help="G-band bright limit", type=float)
  parser.add_argument("-p", action="store_true", dest="pdfOutput", help="Make PDF plot")
  parser.add_argument("-b", action="store_true", dest="pngOutput", help="Make PNG plot")

  args=vars(parser.parse_args())
  return args

if __name__ in ('__main__'):
  args=parseCommandLineArguments()
  plotBrightLimitInV(args['gmag'], pdf=args['pdfOutput'], png=args['pngOutput'])
