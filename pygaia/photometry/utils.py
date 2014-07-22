__all__ = ['vminiFromSpt', 'vabsFromSpt', 'gabsFromSpt']

from scipy import isscalar
from numpy import sqrt, power, amax
from pygaia.photometry.transformations import gminvFromVmini

_sptToVminiVabsDictionary={'B0V':(-0.31, -3.5), 'B1V':(-0.24,-2.7), 'B5V':(-0.08,0.0), 'A0V':(0.01,0.0),
'A5V':(0.16,1.69), 'F0V':(0.38,2.98), 'G0V':(0.67,4.24), 'G2V':(0.72,4.7), 'G5V':(0.74,4.78),
'K0V':(0.87,5.58), 'K1IIIMP':(0.99,1.53), 'K4V':(1.23,7.21), 'K1III':(1.04,2.16), 'M0V':(1.71,8.62),
'M2V':(2.02,9.48), 'M6V':(3.69,14.2), 'M0III':(1.65,-0.66), 'B0I':(-0.22,-7.1), 'B1I':(-0.16,-6.7)}

def vminiFromSpt(spt):
  """
  Obtain (V-I) for the input spectral type.

  Parameters
  ----------
  
  spt - String representing the spectral type of the star.

  Returns
  -------
  
  The value of (V-I).
  """
  if spt in _sptToVminiVabsDictionary:
    return _sptToVminiVabsDictionary[spt][0]
  else:
    message="Unknown spectral type. Allowed values are: "
    for key in _sptToVminiVabsDictionary.keys():
      message += key+" "
    raise Exception(message)

def vabsFromSpt(spt):
  """
  Obtain M_V (absolute magnitude in V-band) for the input spectral type.

  Parameters
  ----------
  
  spt - String representing the spectral type of the star.

  Returns
  -------
  
  The value of M_V.
  """
  if spt in _sptToVminiVabsDictionary:
    return _sptToVminiVabsDictionary[spt][1]
  else:
    message="Unknown spectral type. Allowed values are: "
    for key in _sptToVminiVabsDictionary.keys():
      message += key+" "
    raise Exception(message)

def gabsFromSpt(spt):
  """
  Obtain M_G (absolute magnitude in G-band) for the input spectral type.

  Parameters
  ----------
  
  spt - String representing the spectral type of the star.

  Returns
  -------
  
  The value of M_G.
  """
  if spt in _sptToVminiVabsDictionary:
    return vabsFromSpt(spt) + gminvFromVmini(vminiFromSpt(spt))
  else:
    message="Unknown spectral type. Allowed values are: "
    for key in _sptToVminiVabsDictionary.keys():
      message += key+" "
    raise Exception(message)
