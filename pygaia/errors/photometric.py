__all__ = ['gMagnitudeError', 'bpMagnitudeError', 'rpMagnitudeError']

from numpy import power, sqrt
from utils import calcZ, calcZBpRp

def gMagnitudeError(G, vmini):
  """
  Calculate the single-field-of-view-transit photometric standard errors in the G band as a function
  of G and (V-I).

  Parameters
  ----------

  G     - Value(s) of G-band magnitude.
  vmini - Value(s) of (V-I) colour.

  Returns
  -------

  The G band photometric standard error in units of magnitude.
  """
  z=calcZ(G)
  return 1.0e-3*sqrt(0.02076*z*z + 2.7224*z + 0.004352)

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
  a	=	-0.003201*power(vmini,3) + 0.0589*vmini*vmini + 0.3353*vmini + 0.7927
  b	=	-0.001019*power(vmini,3) + 0.0244*vmini*vmini + 0.1756*vmini + 1.4684
  c	=	-0.004093*power(vmini,3) + 0.0740*vmini*vmini + 0.2834*vmini - 3.4772
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
  a	=	-0.006560*power(vmini,3) + 0.1080*vmini*vmini - 0.6296*vmini + 1.4470
  b	=	-0.003280*power(vmini,3) + 0.0540*vmini*vmini - 0.3148*vmini + 1.7856
  c	=	-0.007992*power(vmini,3) + 0.1482*vmini*vmini - 0.7544*vmini - 3.7232
  return 1.0e-3*sqrt(power(10.0,a)*z*z+power(10.0,b)*z+power(10.0,c))
