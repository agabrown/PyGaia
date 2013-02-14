__all__ = ['enum', 'degreesToRadians', 'radiansToDegrees']

from numpy import pi

def enum(typename, field_names):
  """
  Create a new enumeration type.
  
  Code is copyright (c) Gabriel Genellina, 2010, MIT License.

  Parameters
  ----------

  typename - Name of the enumerated type
  field_names - Names of the fields of the enumerated type
  """

  if isinstance(field_names, str):
    field_names = field_names.replace(',', ' ').split()
  d = dict((reversed(nv) for nv in enumerate(field_names)), __slots__ = ())
  return type(typename, (object,), d)()

def degreesToRadians(angle):
  """
  Convert from degrees to radians.

  Parameters
  ----------

  angle - angle in degrees

  Returns
  -------

  Angle in radians.
  """
  return angle/180.0*pi

def radiansToDegrees(angle):
  """
  Convert from radians to degrees.

  Parameters
  ----------

  angle - angle in radians.

  Returns
  -------

  Angle in degrees.
  """
  return angle/pi*180.0
