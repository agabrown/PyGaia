__all__ = ['CoordinateTransformation', 'Transformations']

from pygaia.astrometry.vectorastrometry import sphericalToCartesian, cartesianToSpherical, \
        elementaryRotationMatrix, normalTriad
from pygaia.utils import enum, degreesToRadians, radiansToDegrees

from numpy import ones_like, array, pi, cos, sin, zeros_like
from numpy import dot, transpose, cross, vstack, diag, sqrt, identity
from numpy.linalg import norm
from scipy import isscalar

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

  def transformProperMotions(self, phi, theta, muphistar, mutheta):
    """
    Converts proper motions from one reference system to another, using the prescriptions in section
    1.5.3 of the Hipparcos Explanatory Volume 1 (equations 1.5.18, 1.5.19).

    Parameters
    ----------

    phi       - The longitude-like angle of the position of the source (radians).
    theta     - The latitude-like angle of the position of the source (radians).
    muphistar - Value of the proper motion in the longitude-like angle, multiplied by cos(latitude).
    mutheta   - Value of the proper motion in the latitude-like angle.

    Returns
    -------

    muphistarrot - Value of the transformed proper motion in the longitude-like angle (including the
    cos(latitude) factor).
    muthetarot   - Value of the transformed proper motion in the latitude-like angle.
    """
    c, s = self._getJacobian(phi,theta)
    return c*muphistar+s*mutheta, c*mutheta-s*muphistar

  def transformSkyCoordinateErrors(self, phi, theta, sigPhiStar, sigTheta, rhoPhiTheta=0):
    """
    Converts the sky coordinate errors from one reference system to another, including the covariance
    term. Equations (1.5.4) and (1.5.20) from section 1.5 in the Hipparcos Explanatory Volume 1 are used.

    Parameters
    ----------

    phi         - The longitude-like angle of the position of the source (radians).
    theta       - The latitude-like angle of the position of the source (radians).
    sigPhiStar  - Standard error in the longitude-like angle of the position of the source (radians or
                  sexagesimal units, including cos(latitude) term)
    sigTheta    - Standard error in the latitude-like angle of the position of the source (radians or
                  sexagesimal units)

    Keywords (optional)
    -------------------

    rhoPhiTheta - Correlation coefficient of the position errors. Set to zero if this keyword is not
                  provided.

    Retuns
    ------
    
    sigPhiRotStar  - The transformed standard error in the longitude-like angle (including
                     cos(latitude) factor)
    sigThetaRot    - The transformed standard error in the latitude-like angle.
    rhoPhiThetaRot - The transformed correlation coefficient.
    """
    if isscalar(rhoPhiTheta) and not isscalar(sigTheta):
      rhoPhiTheta=zeros_like(sigTheta)+rhoPhiTheta
    c, s = self._getJacobian(phi,theta)
    cSqr = c*c
    sSqr = s*s
    covar = sigPhiStar*sigTheta*rhoPhiTheta
    varPhiStar = sigPhiStar*sigPhiStar
    varTheta = sigTheta*sigTheta
    varPhiStarRot = cSqr*varPhiStar+sSqr*varTheta+2.0*covar*c*s
    varThetaRot = sSqr*varPhiStar+cSqr*varTheta-2.0*covar*c*s
    covarRot = (cSqr-sSqr)*covar+c*s*(varTheta-varPhiStar)
    return sqrt(varPhiStarRot), sqrt(varThetaRot), covarRot/sqrt(varPhiStarRot*varThetaRot)

  def transformProperMotionErrors(self, phi, theta, sigMuPhiStar, sigMuTheta, rhoMuPhiMuTheta=0):
    """
    Converts the proper motion errors from one reference system to another, including the covariance
    term. Equations (1.5.4) and (1.5.20) from section 1.5 in the Hipparcos Explanatory Volume 1 are used.

    Parameters
    ----------

    phi             - The longitude-like angle of the position of the source (radians).
    theta           - The latitude-like angle of the position of the source (radians).
    sigMuPhiStar    - Standard error in the proper motion in the longitude-like direction (including
    cos(latitude) factor).
    sigMuTheta      - Standard error in the proper motion in the latitude-like direction.

    Keywords (optional)
    -------------------

    rhoMuPhiMuTheta - Correlation coefficient of the proper motion errors. Set to zero if this 
                      keyword is not provided.

    Retuns
    ------
    
    sigMuPhiRotStar    - The transformed standard error in the proper motion in the longitude direction
    (including cos(latitude) factor).
    sigMuThetaRot      - The transformed standard error in the proper motion in the longitude direction.
    rhoMuPhiMuThetaRot - The transformed correlation coefficient.
    """
    return self.transformSkyCoordinateErrors(phi, theta, sigMuPhiStar, sigMuTheta,
        rhoPhiTheta=rhoMuPhiMuTheta)

  def transformCovarianceMatrix(self, phi, theta, covmat):
      """
      Transform the astrometric covariance matrix to its representation in the new coordinate system.

      Parameters
      ----------

      phi       - The longitude-like angle of the position of the source (radians).
      theta     - The latitude-like angle of the position of the source (radians).
      covmat    - Covariance matrix (5x5) of the astrometric parameters.

      Returns
      -------

      covmat_rot - Covariance matrix in its representation in the new coordinate system.
      """

      c, s = self._getJacobian(phi,theta)
      jacobian = identity(5)
      jacobian[0][0]=c
      jacobian[1][1]=c
      jacobian[3][3]=c
      jacobian[4][4]=c
      jacobian[0][1]=s
      jacobian[1][0]=-s
      jacobian[3][4]=s
      jacobian[4][3]=-s

      return dot( dot(jacobian, covmat), jacobian.transpose() )

  def _getJacobian(self, phi, theta):
    """
    Calculates the Jacobian for the transformation of the position errors and proper motion errors
    between coordinate systems. This Jacobian is also the rotation matrix for the transformation of
    proper motions. See section 1.5.3 of the Hipparcos Explanatory Volume 1 (equation 1.5.20).

    Parameters
    ----------

    phi       - The longitude-like angle of the position of the source (radians).
    theta     - The latitude-like angle of the position of the source (radians).

    Returns
    -------

    jacobian - The Jacobian matrix corresponding to (phi, theta) and the currently desired coordinate
               system transformation.
    """

    p, q, r = normalTriad(phi, theta)

    # zRot = z-axis of new coordinate system expressed in terms of old system
    zRot = self.rotationMatrix[2,:]
    zRotAll = zRot
    if (p.ndim == 2):
      for i in range(p.shape[1]-1):
        zRotAll = vstack((zRotAll,zRot))
    pRot = cross(zRotAll, transpose(r))
    if (p.ndim == 2):
      normPRot = sqrt(diag(dot(pRot,transpose(pRot))))
      for i in range(pRot.shape[0]):
        pRot[i,:] = pRot[i,:]/normPRot[i]
    else:
      pRot = pRot/norm(pRot)

    if (p.ndim == 2):
      return diag(dot(pRot,p)), diag(dot(pRot,q))
    else:
      return dot(pRot,p), dot(pRot,q)
