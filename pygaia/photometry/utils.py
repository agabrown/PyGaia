__all__ = ['vminiFromSpt']

from scipy import isscalar
from numpy import sqrt, power, amax

_sptToVminiDictionary={'B0V':-0.31, 'B1V':-0.24, 'B5V':-0.08, 'A0V':0.01, 'A5V':0.16,
'F0V':0.38, 'G0V':0.67, 'G2V':0.72, 'G5V':0.74, 'K0V':0.87, 'K1IIIMP':0.99, 'K4V':1.23,
'K1III':1.04, 'M6V':3.69}

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
  if spt in _sptToVminiDictionary:
    return _sptToVminiDictionary[spt]
  else:
    message="Unknown spectral type. Allowed values are: "
    for key in _sptToVminiDictionary.keys():
      message += key+" "
    raise Exception(message)
