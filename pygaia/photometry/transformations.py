__all__ = ['gminvFromVmini', 'vminGrvsFromVmini']

from scipy import isscalar
from numpy import sqrt, power, amax

def gminvFromVmini(vmini):
  """
  Calculate the value of (G-V) from (V-I).

  Parameters
  ----------

  vmini - The value of (V-I).

  Returns
  -------

  The value of (G-V)
  """
  return -0.0257-0.0924*vmini-0.1623*vmini*vmini+0.0090*power(vmini,3)

def vminGrvsFromVmini(vmini):
  """
  Calculate (V-Grvs) from (V-I).

  Parameters
  ----------
  
  vmini - The value of (V-I).

  Returns
  -------
  
  The value of (V-Grvs).
  """
  return  0.0119+1.2092*vmini-0.0188*vmini*vmini-0.0005*power(vmini,3)
