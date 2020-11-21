__all__ = ['spherical_to_cartesian', 'cartesian_to_spherical', 'normal_triad', 'elementary_rotation_matrix',
           'phase_space_to_astrometry', 'astrometry_to_phase_space']

import numpy as np

from .constants import au_mas_parsec, au_km_year_per_sec


def spherical_to_cartesian(r, phi, theta):
    """
    Convert spherical to Cartesian coordinates. The input can be scalars or 1-dimensional numpy arrays.
    Note that the angle coordinates follow the astronomical convention of using elevation (declination,
    latitude) rather than its complement (pi/2-elevation), where the latter is commonly used in the
    mathematical treatment of spherical coordinates.

    Parameters
    ----------

    r     - length of input Cartesian vector.
    phi   - longitude-like angle (e.g., right ascension, ecliptic longitude) in radians
    theta - latitide-like angle (e.g., declination, ecliptic latitude) in radians

    Returns
    -------
  
    The Cartesian vector components x, y, z
    """
    ctheta = np.cos(theta)
    x = r * np.cos(phi) * ctheta
    y = r * np.sin(phi) * ctheta
    z = r * np.sin(theta)
    return x, y, z


def cartesian_to_spherical(x, y, z):
    """
    Convert Cartesian to spherical coordinates. The input can be scalars or 1-dimensional numpy arrays.
    Note that the angle coordinates follow the astronomical convention of using elevation (declination,
    latitude) rather than its complement (pi/2-elevation), which is commonly used in the mathematical
    treatment of spherical coordinates.

    Parameters
    ----------
  
    x - Cartesian vector component along the X-axis
    y - Cartesian vector component along the Y-axis
    z - Cartesian vector component along the Z-axis

    Returns
    -------
  
    The spherical coordinates r=np.sqrt(x*x+y*y+z*z), longitude phi, latitude theta.
  
    NOTE THAT THE LONGITUDE ANGLE IS BETWEEN 0 AND +2PI. FOR r=0 AN EXCEPTION IS RAISED.
    """
    rCylSq = x * x + y * y
    r = np.sqrt(rCylSq + z * z)
    if np.any(r == 0.0):
        raise Exception("Error: one or more of the points is at distance zero.")
    phi = np.arctan2(y, x)
    phi = np.where(phi < 0.0, phi + 2 * np.pi, phi)
    return r, phi, np.arctan2(z, np.sqrt(rCylSq))


def normal_triad(phi, theta):
    """
    Calculate the so-called normal triad [p, q, r] which is associated with a spherical coordinate system .
    The three vectors are:

    p - The unit tangent vector in the direction of increasing longitudinal angle phi.
    q - The unit tangent vector in the direction of increasing latitudinal angle theta.
    r - The unit vector toward the point (phi, theta).

    Parameters
    ----------

    phi   - longitude-like angle (e.g., right ascension, ecliptic longitude) in radians
    theta - latitide-like angle (e.g., declination, ecliptic latitude) in radians

    Returns
    -------

    The normal triad as the vectors p, q, r
    """
    sphi = np.sin(phi)
    stheta = np.sin(theta)
    cphi = np.cos(phi)
    ctheta = np.cos(theta)
    p = np.array([-sphi, cphi, np.zeros_like(phi)])
    q = np.array([-stheta * cphi, -stheta * sphi, ctheta])
    r = np.array([ctheta * cphi, ctheta * sphi, stheta])
    return p, q, r


def elementary_rotation_matrix(axis, rotationAngle):
    """
    Construct an elementary rotation matrix describing a rotation around the x, y, or z-axis.

    Parameters
    ----------

    axis          - Axis around which to rotate ("x", "y", or "z")
    rotationAngle - the rotation angle in radians

    Returns
    -------

    The rotation matrix

    Example usage
    -------------

    rotmat = elementaryRotationMatrix("y", pi/6.0)
    """
    if (axis == "x" or axis == "X"):
        return np.array([[1.0, 0.0, 0.0], [0.0, np.cos(rotationAngle), np.sin(rotationAngle)], [0.0,
                                                                                                -np.sin(rotationAngle),
                                                                                                np.cos(rotationAngle)]])
    elif (axis == "y" or axis == "Y"):
        return np.array([[np.cos(rotationAngle), 0.0, -np.sin(rotationAngle)], [0.0, 1.0, 0.0], [np.sin(rotationAngle),
                                                                                                 0.0, np.cos(
                rotationAngle)]])
    elif (axis == "z" or axis == "Z"):
        return np.array([[np.cos(rotationAngle), np.sin(rotationAngle), 0.0], [-np.sin(rotationAngle),
                                                                               np.cos(rotationAngle), 0.0],
                         [0.0, 0.0, 1.0]])
    else:
        raise Exception("Unknown rotation axis " + axis + "!")


def phase_space_to_astrometry(x, y, z, vx, vy, vz):
    """
    From the given phase space coordinates calculate the astrometric observables, including the radial
    velocity, which here is seen as the sixth astrometric parameter. The phase space coordinates are
    assumed to represent barycentric (i.e. centred on the Sun) positions and velocities.

    This function has no mechanism to deal with units. The velocity units are always assumed to be km/s,
    and the code is set up such that for positions in pc, the return units for the astrometry are radians,
    milliarcsec, milliarcsec/year and km/s. For positions in kpc the return units are: radians,
    microarcsec, microarcsec/year, and km/s.

    NOTE that the doppler factor k=1/(1-vrad/c) is NOT used in the calculations. This is not a problem for
    sources moving at typical velocities of Galactic stars.

    Parameters
    ----------

    x - The x component of the barycentric position vector (in pc or kpc).
    y - The y component of the barycentric position vector (in pc or kpc).
    z - The z component of the barycentric position vector (in pc or kpc).
    vx - The x component of the barycentric velocity vector (in km/s).
    vy - The y component of the barycentric velocity vector (in km/s).
    vz - The z component of the barycentric velocity vector (in km/s).

    Returns
    -------

    phi       - The longitude-like angle of the position of the source (radians).
    theta     - The latitude-like angle of the position of the source (radians).
    parallax  - The parallax of the source (in mas or muas, see above)
    muphistar - The proper motion in the longitude-like angle, multiplied by np.cos(theta) (mas/yr or muas/yr,
    see above)
    mutheta   - The proper motion in the latitude-like angle (mas/yr or muas/yr, see above)
    vrad      - The radial velocity (km/s)
    """
    u, phi, theta = cartesian_to_spherical(x, y, z)
    parallax = au_mas_parsec / u
    p, q, r = normal_triad(phi, theta)
    velocitiesArray = np.array([vx, vy, vz])
    if np.isscalar(u):
        muphistar = np.dot(p, velocitiesArray) * parallax / au_km_year_per_sec
        mutheta = np.dot(q, velocitiesArray) * parallax / au_km_year_per_sec
        vrad = np.dot(r, velocitiesArray)
    else:
        muphistar = np.zeros_like(parallax)
        mutheta = np.zeros_like(parallax)
        vrad = np.zeros_like(parallax)
        for i in range(parallax.size):
            muphistar[i] = np.dot(p[:, i], velocitiesArray[:, i]) * parallax[i] / au_km_year_per_sec
            mutheta[i] = np.dot(q[:, i], velocitiesArray[:, i]) * parallax[i] / au_km_year_per_sec
            vrad[i] = np.dot(r[:, i], velocitiesArray[:, i])

    return phi, theta, parallax, muphistar, mutheta, vrad


def astrometry_to_phase_space(phi, theta, parallax, muphistar, mutheta, vrad):
    """
    From the input astrometric parameters calculate the phase space coordinates. The output phase space
    coordinates represent barycentric (i.e. centred on the Sun) positions and velocities.

    This function has no mechanism to deal with units. The code is set up such that for input astrometry
    with parallaxes and proper motions in mas and mas/yr, and radial velocities in km/s, the phase space
    coordinates are in pc and km/s. For input astrometry with parallaxes and proper motions in muas and
    muas/yr, and radial velocities in km/s, the phase space coordinates are in kpc and km/s. Only positive
    parallaxes are accepted, an exception is thrown if this condition is not met.

    NOTE that the doppler factor k=1/(1-vrad/c) is NOT used in the calculations. This is not a problem for
    sources moving at typical velocities of Galactic stars.

    THIS FUNCTION SHOULD NOT BE USED WHEN THE PARALLAXES HAVE RELATIVE ERRORS LARGER THAN ABOUT 20 PER CENT
    (see http://arxiv.org/abs/1507.02105 for example). For astrometric data with relatively large parallax
    errors you should consider doing your analysis in the data space and use forward modelling of some
    kind.

    Parameters
    ----------

    phi       - The longitude-like angle of the position of the source (radians).
    theta     - The latitude-like angle of the position of the source (radians).
    parallax  - The parallax of the source (in mas or muas, see above)
    muphistar - The proper motion in the longitude-like angle, multiplied by np.cos(theta) (mas/yr or muas/yr,
    see above)
    mutheta   - The proper motion in the latitude-like angle (mas/yr or muas/yr, see above)
    vrad      - The radial velocity (km/s)

    Returns
    -------

    x - The x component of the barycentric position vector (in pc or kpc).
    y - The y component of the barycentric position vector (in pc or kpc).
    z - The z component of the barycentric position vector (in pc or kpc).
    vx - The x component of the barycentric velocity vector (in km/s).
    vy - The y component of the barycentric velocity vector (in km/s).
    vz - The z component of the barycentric velocity vector (in km/s).
    """
    if np.any(parallax <= 0.0):
        raise Exception("One or more of the input parallaxes is non-positive")
    x, y, z = spherical_to_cartesian(au_mas_parsec / parallax, phi, theta)
    p, q, r = normal_triad(phi, theta)
    transverseMotionArray = np.array(
        [muphistar * au_km_year_per_sec / parallax, mutheta * au_km_year_per_sec / parallax,
         vrad])
    if np.isscalar(parallax):
        velocityArray = np.dot(np.transpose(np.array([p, q, r])), transverseMotionArray)
        vx = velocityArray[0]
        vy = velocityArray[1]
        vz = velocityArray[2]
    else:
        vx = np.zeros_like(parallax)
        vy = np.zeros_like(parallax)
        vz = np.zeros_like(parallax)
        for i in range(parallax.size):
            velocityArray = np.dot(np.transpose(np.array([p[:, i], q[:, i], r[:, i]])), transverseMotionArray[:, i])
            vx[i] = velocityArray[0]
            vy[i] = velocityArray[1]
            vz[i] = velocityArray[2]
    return x, y, z, vx, vy, vz
