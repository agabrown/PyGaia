__all__ = ['calcZ', 'calcZBpRp', 'calcZAltStartGate', 'averageNumberOfTransits']

from scipy import isscalar
from numpy import sqrt, power, amax, array, floor, sin
from numpy import genfromtxt
from pkg_resources import resource_stream

_table = resource_stream('pygaia', 'data/errorFactorVariationBeta.txt')
_averageTransitNumber = genfromtxt(_table,
    skip_header=16, skip_footer=1,
    names=['sinBeta','betaMin','betaMax','nTransits',
      'alphaStar','delta','parallax','muAlphaStar','muDelta'])['nTransits']
_numStepsSinBeta = len(_averageTransitNumber)

def calcZ(G):
  """
  Calculate the value for the parameter z in the formula for parallax and G magnitude errors as a
  function of G and (V-I).

  Parameters
  ----------
  
  G - Value of G-band magnitude.

  Returns
  -------
  
  Value of z.
  """
  gatefloor=power(10.0,0.4*(12.0-15.0))
  if isscalar(G):
   result=amax((gatefloor,power(10.0,0.4*(G-15.0))))
  else :
    result=power(10.0,0.4*(G-15.0))
    indices=(result<gatefloor)
    result[indices]=gatefloor
  return result

def calcZBpRp(G):
  """
  Calculate the value for the parameter z in the formula for the BP and RP magnitude errors as a
  function of G and (V-I).

  Parameters
  ----------
  
  G - Value of G-band magnitude.

  Returns
  -------
  
  Value of z for BP/RP.
  """
  gatefloor=power(10.0,0.4*(11.0-15.0))
  if isscalar(G):
   result=amax((gatefloor,power(10.0,0.4*(G-15.0))))
  else :
    result=power(10.0,0.4*(G-15.0))
    indices=(result<gatefloor)
    result[indices]=gatefloor
  return result

def calcZAltStartGate(G):
  """
  Calculate the value of z in the formula for the parallax errors. In this case assume gating starts at
  G=13.3 (to simulate bright star worst performance)

  Parameters
  ----------
  
  G - Value of G-band magnitude.

  Returns
  -------
  
  Value of z.
  """
  gatefloor=power(10.0,0.4*(13.3-15.0))
  if isscalar(G):
   result=amax((gatefloor,power(10.0,0.4*(G-15.0))))
  else :
    result=power(10.0,0.4*(G-15.0))
    indices=(result<gatefloor)
    result[indices]=gatefloor
  return result

def averageNumberOfTransits(beta):
  """
  Returns the number of transits across the Gaia focal plane averaged over ecliptic longitude.

  Parameters
  ----------

  beta - Value(s) of the Ecliptic latitude.

  Returns
  -------

  Average number of transits for the input values of beta.
  """
  indices = array(floor(abs(sin(beta))*_numStepsSinBeta), dtype=int)
  indices[(indices==_numStepsSinBeta)] = _numStepsSinBeta-1
  return _averageTransitNumber[indices]
