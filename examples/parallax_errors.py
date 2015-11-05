#!/usr/bin/env python

# Predict parallax errors as a function of G and (V-I) according to the
# interpolation formula on the Gaia science performance web-pages.

import numpy as np
from pygaia.photometry.transformations import gminvFromVmini
from pygaia.errors.astrometric import parallaxErrorSkyAvg

from os import environ as env
try:
  import argparse
except ImportError:
  raise ImportError("""The argparse module does not seem to be available.\n Either upgrade to python 2.7
  or adapt the script to using the optparse module instead.""")

def calcParallaxError(args):
  """
  Calculate the parallax error for the given input source magnitude and colour.

  :argument args: command line arguments
  """
  gmag=float(args['gmag'])
  vmini=float(args['vmini'])
  sigmaPar=parallaxErrorSkyAvg(gmag, vmini)
  gminv=gminvFromVmini(vmini)
  print("G = {0}".format(gmag))
  print("V = {0}".format(gmag-gminv))
  print("(V-I) = {0}".format(vmini))
  print("(G-V) = {0}".format(gminv))
  print("standard error = {0} muas".format(sigmaPar))

def parseCommandLineArguments():
  """
  Set up command line parsing.
  """
  parser = argparse.ArgumentParser(description="Calculate parallax error for given G and (V-I)")
  parser.add_argument("gmag", help="G-band magnitude of source", type=float)
  parser.add_argument("vmini", help="(V-I) colour of source", type=float)

  args=vars(parser.parse_args())
  return args

if __name__ in ('__main__'):
  args=parseCommandLineArguments()
  calcParallaxError(args)
