__all__ = ['gMagnitudeError', 'bpMagnitudeError', 'rpMagnitudeError']

from numpy import power, sqrt
from utils import calcZ, calcZBpRp

def gMagnitudeError(G):
  """
  Calculate the single-field-of-view-transit photometric standard errors in the G band as a function
  of G.

  Parameters
  ----------

  G     - Value(s) of G-band magnitude.

  Returns
  -------

  The G band photometric standard error in units of magnitude.
  """
  z=calcZ(G)
  return 1.0e-3*sqrt(0.04895*z*z + 1.8633*z + 0.0001985)

def bpMagnitudeError(G, vmini):
  """
  Calculate the single-field-of-view-transit photometric standard errors in the BP band as a function
  of G and (V-I). Note: this refers to the integrated flux from the BP spectrophotometer.

  Parameters
  ----------

  G     - Value(s) of G-band magnitude.
  vmini - Value(s) of (V-I) colour.

  Returns
  -------

  The BP band photometric standard error in units of magnitude.
  """
  z=calcZBpRp(G)
  a	=	-0.000562*power(vmini,3) + 0.044390*vmini*vmini + 0.355123*vmini + 1.043270
  b	=	-0.000400*power(vmini,3) + 0.018878*vmini*vmini + 0.195768*vmini + 1.465592
  c	=	+0.000262*power(vmini,3) + 0.060769*vmini*vmini - 0.205807*vmini - 1.866968
  return 1.0e-3*sqrt(power(10.0,a)*z*z+power(10.0,b)*z+power(10.0,c))

def rpMagnitudeError(G, vmini):
  """
  Calculate the single-field-of-view-transit photometric standard errors in the RP band as a function
  of G and (V-I). Note: this refers to the integrated flux from the RP spectrophotometer.

  Parameters
  ----------

  G     - Value(s) of G-band magnitude.
  vmini - Value(s) of (V-I) colour.

  Returns
  -------

  The RP band photometric standard error in units of magnitude.
  """
  z=calcZBpRp(G)
  a	=	-0.007597*power(vmini,3) + 0.114126*vmini*vmini - 0.636628*vmini + 1.615927
  b	=	-0.003803*power(vmini,3) + 0.057112*vmini*vmini - 0.318499*vmini + 1.783906
  c	=	-0.001923*power(vmini,3) + 0.027352*vmini*vmini - 0.091569*vmini - 3.042268
  return 1.0e-3*sqrt(power(10.0,a)*z*z+power(10.0,b)*z+power(10.0,c))
