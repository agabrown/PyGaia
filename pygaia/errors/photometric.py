__all__ = ['gMagnitudeError', 'bpMagnitudeError', 'rpMagnitudeError']

from numpy import power, sqrt
from pygaia.errors.utils import calcZ, calcZBpRp

#
# Margin to include on predicted standard errors (i.e. multiply prediction by this value).
_scienceMargin = 1.2

#
# Mean number of CCDs crossed by a source in the AF field (G-band photometry)
_meanNumCcds = (7.0*9.0-1.0)/7.0

#
# End-of-mission CCD transit calibration floor on the photometric errors
_eomCalibrationFloorG = 3.0e-3/sqrt(_meanNumCcds)
_eomCalibrationFloorBP = 5.0e-3
_eomCalibrationFloorRP = 5.0e-3

def gMagnitudeError(G):
  """
  Calculate the single-field-of-view-transit photometric standard error in the G band as a function
  of G. A 20% margin is included.

  Parameters
  ----------

  G     - Value(s) of G-band magnitude.

  Returns
  -------

  The G band photometric standard error in units of magnitude.
  """
  z=calcZ(G)
  return 1.0e-3*sqrt(0.04895*z*z + 1.8633*z + 0.0001985) * _scienceMargin

def gMagnitudeErrorEoM(G, nobs=70):
  """
  Calculate the end of mission photometric standard error in the G band as a function
  of G. A 20% margin is included.

  Parameters
  ----------

  G     - Value(s) of G-band magnitude.

  Keywords
  --------

  nobs  - Number of observations collected (default 70).

  Returns
  -------

  The G band photometric standard error in units of magnitude.
  """
  return sqrt( (power(gMagnitudeError(G)/_scienceMargin,2) +
      _eomCalibrationFloorG*_eomCalibrationFloorG)/nobs ) * _scienceMargin

def bpMagnitudeError(G, vmini):
  """
  Calculate the single-field-of-view-transit photometric standard error in the BP band as a function
  of G and (V-I). Note: this refers to the integrated flux from the BP spectrophotometer. A margin of 20%
  is included.

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

def bpMagnitudeErrorEoM(G, vmini, nobs=70):
  """
  Calculate the end-of-mission photometric standard error in the BP band as a function of G and (V-I).
  Note: this refers to the integrated flux from the BP spectrophotometer. A margin of 20% is included.

  Parameters
  ----------

  G     - Value(s) of G-band magnitude.
  vmini - Value(s) of (V-I) colour.

  Keywords
  --------

  nobs  - Number of observations collected (default 70).

  Returns
  -------

  The BP band photometric standard error in units of magnitude.
  """
  return sqrt( (power(bpMagnitudeError(G, vmini)/_scienceMargin,2) +
      _eomCalibrationFloorBP*_eomCalibrationFloorBP)/nobs ) * _scienceMargin

def rpMagnitudeError(G, vmini):
  """
  Calculate the single-field-of-view-transit photometric standard error in the RP band as a function
  of G and (V-I). Note: this refers to the integrated flux from the RP spectrophotometer. A margin of 20%
  is included.

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

def rpMagnitudeErrorEoM(G, vmini, nobs=70):
  """
  Calculate the end-of-mission photometric standard error in the RP band as a function of G and (V-I).
  Note: this refers to the integrated flux from the RP spectrophotometer. A margin of 20% is included.

  Parameters
  ----------

  G     - Value(s) of G-band magnitude.
  vmini - Value(s) of (V-I) colour.

  Keywords
  --------

  nobs  - Number of observations collected (default 70).

  Returns
  -------

  The RP band photometric standard error in units of magnitude.
  """
  return sqrt( (power(rpMagnitudeError(G, vmini)/_scienceMargin,2) +
      _eomCalibrationFloorRP*_eomCalibrationFloorRP)/nobs ) * _scienceMargin

