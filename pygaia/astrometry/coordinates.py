__all__ = ['CoordinateTransformation', 'Transformations']

from vectorastrometry import sphericalToCartesian, cartesianToSpherical, elementaryRotationMatrix
from pygaia.utils import enum, degreesToRadians, radiansToDegrees

from numpy import ones_like, array, pi, cos, sin
from numpy import dot, transpose

# Obliquity of the Ecliptic (arcsec)
_obliquityOfEcliptic = degreesToRadians(84381.41100/3600.0)

# Galactic pole in ICRS coordinates (see Hipparcos Explanatory Vol 1 section 1.5, and Murray, 1983,
# section 10.2)
_alphaGalPole = degreesToRadians(192.85948)
_deltaGalPole = degreesToRadians(27.12825)
# The galactic longitude of the ascending node of the galactic plane on the equator of ICRS (see
# Hipparcos Explanatory Vol 1 section 1.5, and Murray, 1983, section 10.2)
_omega = degreesToRadians(32.93192)

# Rotation matrix for the transformation from ICRS to Galactic coordinates. See equation (4.25) in
# chapter 4.5 of "Astrometry for Astrophysics", 2012, van Altena et al.
_matA = elementaryRotationMatrix("z",pi/2.0+_alphaGalPole)
_matB = elementaryRotationMatrix("x",pi/2.0-_deltaGalPole)
_matC = elementaryRotationMatrix("z",-_omega)
_rotationMatrixIcrsToGalactic=dot(_matC,dot(_matB,_matA))

# Alternative way to calculate the rotation matrix from ICRS to Galactic coordinates. First calculate
# the vectors describing the Galactic coordinate reference frame expressed within the ICRS.
# _vecN = array([0,0,1])
# _vecG3 = array([cos(_alphaGalPole)*cos(_deltaGalPole), sin(_alphaGalPole)*cos(_deltaGalPole),
#   sin(_deltaGalPole)])
# _vecG0 = cross(_vecN,_vecG3)
# _vecG0 = _vecG0/sqrt(dot(_vecG0,_vecG0))
# _vecG1 = -sin(_omega)*cross(_vecG3,_vecG0)+cos(_omega)*_vecG0
# _vecG2 = cross(_vecG3,_vecG1)
# _rotationMatrixIcrsToGalactic=array([_vecG1,_vecG2,_vecG3])

# Rotation matrix for the transformation from Galactic to ICRS coordinates.
_rotationMatrixGalacticToIcrs = transpose(_rotationMatrixIcrsToGalactic)

# Rotation matrix for the transformation from Ecliptic to ICRS coordinates.
_rotationMatrixEclipticToIcrs = elementaryRotationMatrix("x", -1*(_obliquityOfEcliptic))

# Rotation matrix for the transformation from ICRS to Ecliptic coordinates.
_rotationMatrixIcrsToEcliptic = transpose(_rotationMatrixEclipticToIcrs)

# Rotation matrix for the transformation from Galactic to Ecliptic coordinates.
_rotationMatrixGalacticToEcliptic = dot(_rotationMatrixIcrsToEcliptic,_rotationMatrixGalacticToIcrs)

# Rotation matrix for the transformation from Ecliptic to Galactic coordinates.
_rotationMatrixEclipticToGalactic = transpose(_rotationMatrixGalacticToEcliptic)

Transformations = enum('Transformations', ['GAL2ICRS', 'ICRS2GAL', 'ECL2ICRS', 'ICRS2ECL', 'GAL2ECL', 'ECL2GAL'])

_rotationMatrixMap = {Transformations.GAL2ICRS:_rotationMatrixGalacticToIcrs,
    Transformations.ICRS2GAL:_rotationMatrixIcrsToGalactic,
    Transformations.ECL2ICRS:_rotationMatrixEclipticToIcrs,
    Transformations.ICRS2ECL:_rotationMatrixIcrsToEcliptic,
    Transformations.GAL2ECL:_rotationMatrixGalacticToEcliptic,
    Transformations.ECL2GAL:_rotationMatrixEclipticToGalactic}

_transformationStringMap = {Transformations.GAL2ICRS:("galactic", "ICRS"),
    Transformations.ICRS2GAL:("ICRS","galactic"),
    Transformations.ECL2ICRS:("ecliptic","ICRS"),
    Transformations.ICRS2ECL:("ICRS","ecliptic"),
    Transformations.GAL2ECL:("galactic","ecliptic"),
    Transformations.ECL2GAL:("ecliptic","galactic")}

class CoordinateTransformation:
  """
  Provides methods for carrying out transformations between different coordinate (reference) systems.
  Currently the following are supported:

  ICRS or Equatorial (right ascension, declination)
  Galactic (longitude, latitude)
  Ecliptic (longitude, latitude)

  The transformations can be applied to sky coordinates or Cartesian coordinates.
  """

  def __init__(self, desiredTransformation):
    """
    Class constructor/initializer. Sets the type of transformation to be used.

    Parameters
    ----------

    desiredTransformation - The kind of coordinate transformation that should be provided. For example
    Transformations.GAL2ECL
    """
    self.rotationMatrix=_rotationMatrixMap[desiredTransformation]
    self.transformationStrings=_transformationStringMap[desiredTransformation]

  def transformCartesianCoordinates(self, x, y, z):
    """
    Rotates Cartesian coordinates from one reference system to another using the rotation matrix with
    which the class was initialized. The inputs  can be scalars or 1-dimensional numpy arrays.

    Parameters
    ----------

    x - Value of X-coordinate in original reference system
    y - Value of Y-coordinate in original reference system
    z - Value of Z-coordinate in original reference system

    Returns
    -------

    xrot - Value of X-coordinate after rotation
    yrot - Value of Y-coordinate after rotation
    zrot - Value of Z-coordinate after rotation
    """
    xrot, yrot, zrot = dot(self.rotationMatrix,[x,y,z])
    return xrot, yrot, zrot

  def transformSkyCoordinates(self, phi, theta):
    """
    Converts sky coordinates from one reference system to another, making use of the rotation matrix with
    which the class was initialized. Inputs can be scalars or 1-dimensional numpy arrays.

    Parameters
    ----------

    phi   - Value of the azimuthal angle (right ascension, longitude) in radians.
    theta - Value of the elevation angle (declination, latitude) in radians.

    Returns
    -------

    phirot   - Value of the transformed azimuthal angle in radians.
    thetarot - Value of the transformed elevation angle in radians.
    """
    r=ones_like(phi)
    x, y, z = sphericalToCartesian(r, phi, theta)
    xrot, yrot, zrot = self.transformCartesianCoordinates(x, y, z)
    r, phirot, thetarot = cartesianToSpherical(xrot, yrot, zrot)
    return phirot, thetarot
