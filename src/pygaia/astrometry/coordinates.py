"""
Provides tools for coordinate transformation and epoch propagation.
"""
import numpy as np
from enum import Enum, auto

from pygaia.astrometry.constants import au_km_year_per_sec
from pygaia.astrometry.vectorastrometry import (
    spherical_to_cartesian,
    cartesian_to_spherical,
    elementary_rotation_matrix,
    normal_triad,
)

__all__ = [
    "CoordinateTransformation",
    "Transformations",
    "EpochPropagation",
    "angular_distance",
]

# Obliquity of the Ecliptic (arcsec expressed in radians)
_obliquityOfEcliptic = np.deg2rad(84381.41100 / 3600.0)

# Galactic pole in ICRS coordinates (see Hipparcos Explanatory Vol 1 section 1.5, and
# Murray, 1983, # section 10.2)
_alphaGalPole = np.deg2rad(192.85948)
_deltaGalPole = np.deg2rad(27.12825)

# The galactic longitude of the ascending node of the galactic plane on the equator of
# ICRS (see Hipparcos Explanatory Vol 1 section 1.5, and Murray, 1983, section 10.2)
_omega = np.deg2rad(32.93192)

# Rotation matrix for the transformation from ICRS to Galactic coordinates. See equation
# (4.25) in chapter 4.5 of "Astrometry for Astrophysics", 2012, van Altena et al.
_matA = elementary_rotation_matrix("z", np.pi / 2.0 + _alphaGalPole)
_matB = elementary_rotation_matrix("x", np.pi / 2.0 - _deltaGalPole)
_matC = elementary_rotation_matrix("z", -_omega)
_rotationMatrixIcrsToGalactic = np.dot(_matC, np.dot(_matB, _matA))

# Alternative way to calculate the rotation matrix from ICRS to Galactic coordinates.
# First calculate the vectors describing the Galactic coordinate reference frame
# expressed within the ICRS.
#
# _vecN = array([0,0,1])
# _vecG3 = array([np.cos(_alphaGalPole)*np.cos(_deltaGalPole),
#       np.sin(_alphaGalPole)*np.cos(_deltaGalPole), np.sin(_deltaGalPole)])
# _vecG0 = np.cross(_vecN,_vecG3)
# _vecG0 = _vecG0/np.sqrt(np.dot(_vecG0,_vecG0))
# _vecG1 = -np.sin(_omega)*np.cross(_vecG3,_vecG0)+np.cos(_omega)*_vecG0
# _vecG2 = np.cross(_vecG3,_vecG1)
# _rotationMatrixIcrsToGalactic=array([_vecG1,_vecG2,_vecG3])

# Rotation matrix for the transformation from Galactic to ICRS coordinates.
_rotationMatrixGalacticToIcrs = np.transpose(_rotationMatrixIcrsToGalactic)

# Rotation matrix for the transformation from Ecliptic to ICRS coordinates.
_rotationMatrixEclipticToIcrs = elementary_rotation_matrix(
    "x", -1 * _obliquityOfEcliptic
)

# Rotation matrix for the transformation from ICRS to Ecliptic coordinates.
_rotationMatrixIcrsToEcliptic = np.transpose(_rotationMatrixEclipticToIcrs)

# Rotation matrix for the transformation from Galactic to Ecliptic coordinates.
_rotationMatrixGalacticToEcliptic = np.dot(
    _rotationMatrixIcrsToEcliptic, _rotationMatrixGalacticToIcrs
)

# Rotation matrix for the transformation from Ecliptic to Galactic coordinates.
_rotationMatrixEclipticToGalactic = np.transpose(_rotationMatrixGalacticToEcliptic)


class Transformations(Enum):
    """
    Enumeration with the available coordinate tranformations.
    """

    GAL2ICRS = ("Galactic", "ICRS")
    """
    Transform from Galactic to ICRS coordinates.
    """
    ICRS2GAL = ("ICRS", "Galactic")
    """
    Transform from ICRS to Galactic coordinates.
    """
    ECL2ICRS = ("Ecliptic", "ICRS")
    """
    Transform from Ecliptic to ICRS coordinates.
    """
    ICRS2ECL = ("ICRS", "Ecliptic")
    """
    Transform from ICRS to Ecliptic coordinates.
    """
    GAL2ECL = ("Galactic", "Ecliptic")
    """
    Transform from Galactic to Ecliptic coordinates.
    """
    ECL2GAL = ("Ecliptic", "Galactic")
    """
    Transform from Ecliptic to Galactic coordinates.
    """

    def __init__(self, fromsys, tosys):
        self.fromsystem = fromsys
        self.tosystem = tosys

    @property
    def from_system(self):
        """
        Coordinate system from which the transformation starts.

        Returns
        -------
        fromsystem : str
            String representation of the starting system.
        """
        return self.fromsystem

    @property
    def to_system(self):
        """
        Target coordinate system of the transformation.

        Returns
        -------
        tosystem : str
            String representation of the target system.
        """
        return self.tosystem


_rotation_matrix_dict = {
    Transformations.GAL2ICRS: _rotationMatrixGalacticToIcrs,
    Transformations.ICRS2GAL: _rotationMatrixIcrsToGalactic,
    Transformations.ECL2ICRS: _rotationMatrixEclipticToIcrs,
    Transformations.ICRS2ECL: _rotationMatrixIcrsToEcliptic,
    Transformations.GAL2ECL: _rotationMatrixGalacticToEcliptic,
    Transformations.ECL2GAL: _rotationMatrixEclipticToGalactic,
}


def angular_distance(phi1, theta1, phi2, theta2, return_posangle=False):
    """
    Calculate the angular distance between pairs of sky coordinates.

    Parameters
    ----------
    phi1 : float
        Longitude of first coordinate (radians).
    theta1 : float
        Latitude of first coordinate (radians).
    phi2 : float
        Longitude of second coordinate (radians).
    theta2 : float
        Latitude of second coordinate (radians).
    return_posangle : boolean
        If true return the position angle of the point (phi2, theta2), as defined in the
        tangent plane at (phi1, theta1). Where the position angle is measured from local
        north over east and ranges between 0 and 2*pi. If the distance between the two
        points is 0 or pi to within floating point precision (as given by
        np.finfo(float).resolution) then 0 is returned for the position angle. This is a
        choice for the cases where mathematically the position angle is not defined.
        This resolution (1e-15) corresponds to about 0.2 nanoarcseconds.

    Returns
    -------
    angdist : float
        Angular distance in radians.
    posangle : float
        Position angle as defined above in radians.
    """
    # The Formula below is numerically more stable than np.arccos( np.sin(theta1)*np.sin(theta2) + np.cos(phi2-phi1)*np.cos(theta1)*np.cos(theta2) ). See: https://en.wikipedia.org/wiki/Great-circle_distance
    dist = np.arctan2(
        np.sqrt(
            (np.cos(theta2) * np.sin(phi2 - phi1)) ** 2
            + (
                np.cos(theta1) * np.sin(theta2)
                - np.sin(theta1) * np.cos(theta2) * np.cos(phi2 - phi1)
            )
            ** 2
        ),
        np.sin(theta1) * np.sin(theta2)
        + np.cos(phi2 - phi1) * np.cos(theta1) * np.cos(theta2),
    )
    if not return_posangle:
        return dist
    else:
        posangle = np.arctan2(
            np.cos(theta2) * np.sin(phi2 - phi1),
            np.cos(theta1) * np.sin(theta2)
            - np.sin(theta1) * np.cos(theta2) * np.cos(phi2 - phi1),
        )
        posangle = np.where(
            np.isclose(np.sin(dist), 0, atol=np.finfo(float).resolution), 0.0, posangle
        )
        return dist, np.squeeze(np.where(posangle < 0, posangle + 2 * np.pi, posangle))


class CoordinateTransformation:
    """
    Provides methods for carrying out transformations between different coordinate
    (reference) systems.  Transformations between the following systems are supported:

    * ICRS or Equatorial (right ascension, declination)
    * Galactic (longitude, latitude)
    * Ecliptic (longitude, latitude)

    The transformations can be applied to sky coordinates or Cartesian coordinates.

    .. note::
        The transformations here are only intended for systems that represent a rotated
        version of the ICRS. For transformations from and to non-ICRS systems (such as
        FK5) refer to Astropy's `coordinates
        <https://docs.astropy.org/en/stable/coordinates/index.html>`_ package.

    Attributes
    ----------
    transformation : Transformations
        The transformation represented by the instance of this class.
    rotationMatrix : array
        Rotation matrix corresponding to the transformation.
    """

    def __init__(self, desired_transformation):
        """
        Class constructor/initializer. Sets the type of transformation to be used.

        Parameters
        ----------
        desiredTransformation : Transformations
            The kind of coordinate transformation that should be provided. For example
            Transformations.GAL2ECL.
        """
        self.transformation = desired_transformation
        self.rotationMatrix = _rotation_matrix_dict[desired_transformation]

    def start_coordinates(self):
        """
        Query which is the starting coordinate system of this transformation.

        Returns
        -------
        from_sys : str
            String with starting coordinate system.
        """
        return self.transformation.from_system

    def target_coordinates(self):
        """
        Query which is the target coordinate system of this transformation.

        Returns
        -------
        to_sys : str
            String with target coordinate system.
        """
        return self.transformation.to_system

    def transform_cartesian_coordinates(self, x, y, z):
        """
        Rotate Cartesian coordinates from one reference system to another.

        Use the rotation matrix with which the class was initialized. The inputs can be
        scalars or 1-dimensional numpy arrays. The length scale units of the inputs can
        be arbitrary but should of course be consistent.

        Parameters
        ----------
        x : float or array
            Value of X-coordinate in original reference system.
        y : float or array
            Value of Y-coordinate in original reference system.
        z : float or array
            Value of Z-coordinate in original reference system.

        Returns
        -------
        xrot : float or array
            Value of X-coordinate after rotation
        yrot : float or array
            Value of Y-coordinate after rotation
        zrot : float or array
            Value of Z-coordinate after rotation
        """
        xrot, yrot, zrot = np.dot(self.rotationMatrix, [x, y, z])
        return xrot, yrot, zrot

    def transform_sky_coordinates(self, phi, theta):
        """
        Convert sky coordinates from one reference system to another.

        Use of the rotation matrix with which the class was initialized. Inputs can be
        scalars or 1-dimensional numpy arrays.

        Parameters
        ----------
        phi : float or array
            Value of the longitude-like angle (right ascension, longitude) in radians.
        theta : float or array
            Value of the latitude-like angle (declination, latitude) in radians.

        Returns
        -------
        phirot : float or array
            Value of the transformed longitude-like angle in radians.
        thetarot : float or array
            Value of the transformed latitude-like angle in radians.
        """
        r = np.ones_like(phi)
        x, y, z = spherical_to_cartesian(r, phi, theta)
        xrot, yrot, zrot = self.transform_cartesian_coordinates(x, y, z)
        r, phirot, thetarot = cartesian_to_spherical(xrot, yrot, zrot)
        return phirot, thetarot

    def transform_proper_motions(self, phi, theta, muphistar, mutheta):
        """
        Converts proper motions from one reference system to another.

        Use the prescriptions in section 1.5.3 of the `Hipparcos Explanatory Volume 1
        <https://www.cosmos.esa.int/documents/532822/552851/vol1_all.pdf/99adf6e3-6893-4824-8fc2-8d3c9cbba2b5>`_
        (equations 1.5.18, 1.5.19).

        Parameters
        ----------
        phi : float or array
            The longitude-like angle of the position of the source (radians).
        theta : float or array
            The latitude-like angle of the position of the source (radians).
        muphistar : float or array
            Value of the proper motion in the longitude-like angle, multiplied by
            cos(latitude).
        mutheta : float or array
            Value of the proper motion in the latitude-like angle.

        Returns
        -------
        muphistarrot : float or array
            Value of the transformed proper motion in the longitude-like angle
            (including the cos(latitude) factor).
        muthetarot : float or array
            Value of the transformed proper motion in the latitude-like angle.
        """
        c, s = self._get_jacobian(phi, theta)
        return c * muphistar + s * mutheta, c * mutheta - s * muphistar

    def transform_sky_coordinate_errors(
        self, phi, theta, sigphistar, sigtheta, rho_phi_theta=0
    ):
        """
        Converts the sky coordinate uncertainties from one reference system to another,
        including the covariance term.

        Equations (1.5.4) and (1.5.20) from section 1.5 in the `Hipparcos Explanatory
        Volume 1
        <https://www.cosmos.esa.int/documents/532822/552851/vol1_all.pdf/99adf6e3-6893-4824-8fc2-8d3c9cbba2b5>`_
        are used.

        Parameters
        ----------
        phi : float or array
            The longitude-like angle of the position of the source (radians).
        theta : float or array
            The latitude-like angle of the position of the source (radians).
        sigphistar : float or array
            Uncertainty in the longitude-like angle of the position of the source
            (radians or sexagesimal units, including cos(latitude) term)
        sigtheta : float or array
            Uncertainty in the latitude-like angle of the position of the source
            (radians or sexagesimal units)
        rho_phi_theta : float or array
            Correlation coefficient of the position uncertainties. Zero by default

        Returns
        -------
        sigPhiRotStar : float or array
            The transformed uncertainty in the longitude-like angle (including
            cos(latitude) factor)
        sigThetaRot : float or array
            The transformed uncertainty in the latitude-like angle.
        rhoPhiThetaRot : float or array
            The transformed correlation coefficient.
        """
        if np.isscalar(rho_phi_theta) and not np.isscalar(sigtheta):
            rho_phi_theta = np.zeros_like(sigtheta) + rho_phi_theta
        c, s = self._get_jacobian(phi, theta)
        csqr = c * c
        ssqr = s * s
        covar = sigphistar * sigtheta * rho_phi_theta
        varphistar = sigphistar * sigphistar
        vartheta = sigtheta * sigtheta
        varphistar_rot = csqr * varphistar + ssqr * vartheta + 2.0 * covar * c * s
        vartheta_rot = ssqr * varphistar + csqr * vartheta - 2.0 * covar * c * s
        covar_rot = (csqr - ssqr) * covar + c * s * (vartheta - varphistar)
        return (
            np.sqrt(varphistar_rot),
            np.sqrt(vartheta_rot),
            covar_rot / np.sqrt(varphistar_rot * vartheta_rot),
        )

    def transform_proper_motion_errors(
        self, phi, theta, sigmuphistar, sigmutheta, rho_muphi_mutheta=0
    ):
        """
        Converts the proper motion uncertainties from one reference system to another,
        including the covariance term.

        Equations (1.5.4) and (1.5.20) from section 1.5 in the `Hipparcos Explanatory
        Volume 1
        <https://www.cosmos.esa.int/documents/532822/552851/vol1_all.pdf/99adf6e3-6893-4824-8fc2-8d3c9cbba2b5>`_
        are used.

        Parameters
        ----------
        phi : float or array
            The longitude-like angle of the position of the source (radians).
        theta : float or array
            The latitude-like angle of the position of the source (radians).
        sigmuphistar : float or attay
            Uncertainty in the proper motion in the longitude-like direction
            (including cos(latitude) factor).
        sigmutheta : float or array
            Uncertainty in the proper motion in the latitude-like direction.
        rho_muphi_mutheta : float or array
            Correlation coefficient of the proper motion errors. Zero by default.

        Returns
        -------
        sigMuPhiRotStar : float or array
            The transformed uncertainty in the proper motion in the longitude direction
            (including np.cos(latitude) factor).
        sigMuThetaRot : float or array
            The transformed uncertainty in the proper motion in the longitude direction.
        rhoMuPhiMuThetaRot : float or array
            The transformed correlation coefficient.
        """
        return self.transform_sky_coordinate_errors(
            phi, theta, sigmuphistar, sigmutheta, rho_phi_theta=rho_muphi_mutheta
        )

    def transform_covariance_matrix(self, phi, theta, covmat):
        """
        Transform the astrometric covariance matrix to its representation in the new
        coordinate system.

        Parameters
        ----------
        phi : float
            The longitude-like angle of the position of the source (radians).
        theta : float
            The latitude-like angle of the position of the source (radians).
        covmat : array
            Covariance matrix (5,5) of the astrometric parameters.

        Returns
        -------
        covmat_rot : array
            Covariance matrix in its representation in the new coordinate system.
        """

        c, s = self._get_jacobian(phi, theta)
        jacobian = np.identity(5)
        jacobian[0][0] = c
        jacobian[1][1] = c
        jacobian[3][3] = c
        jacobian[4][4] = c
        jacobian[0][1] = s
        jacobian[1][0] = -s
        jacobian[3][4] = s
        jacobian[4][3] = -s

        return np.dot(np.dot(jacobian, covmat), jacobian.T)

    def _get_jacobian(self, phi, theta):
        r"""
        Calculates the Jacobian for the transformation of the position and proper motion
        uncertainties between coordinate systems.

        This Jacobian is also the rotation matrix for the transformation of proper
        motions. See section 1.5.3 of the `Hipparcos Explanatory Volume 1
        <https://www.cosmos.esa.int/documents/532822/552851/vol1_all.pdf/99adf6e3-6893-4824-8fc2-8d3c9cbba2b5>`_
        (equation 1.5.20). This matrix has the following form:

        .. math::

            J = \begin{pmatrix}
                c & s \\
                -s & c \\
            \end{pmatrix}

        Parameters
        ----------
        phi : float or array
            The longitude-like angle of the position of the source (radians).
        theta : float or array
            The latitude-like angle of the position of the source (radians).

        Returns
        -------
        c, s : float or array
            The Jacobian matrix elements c and s corresponding to (phi, theta) and the
            currently desired coordinate system transformation.
        """

        p, q, r = normal_triad(phi, theta)

        # z_rot = z-axis of new coordinate system expressed in terms of old system
        z_rot = self.rotationMatrix[2, :]

        if p.ndim == 2:
            N = p.shape[1]
            z_rot_all = np.tile(z_rot, N).reshape(N, 3)
            p_rot = np.cross(z_rot_all, r.T)
            p_rot = p_rot / np.linalg.norm(p_rot, axis=1)[:, None]
            c = np.sum(p_rot * p.T, axis=1)
            s = np.sum(p_rot * q.T, axis=1)

            return c, s
        else:
            p_rot = np.cross(z_rot, r.T)
            p_rot = p_rot / np.linalg.norm(p_rot)
            return np.dot(p_rot, p), np.dot(p_rot, q)


class EpochPropagation:
    """
    Provides methods for transforming the astrometry and radial velocity of a given
    source to a different epoch.

    The formulae for rigorous epoch transformation can be found on the `Gaia
    documentation pages
    <https://gea.esac.esa.int/archive/documentation/GDR3/Data_processing/chap_cu3ast/sec_cu3ast_intro/ssec_cu3ast_intro_tansforms.html>`_.

    Attributes
    ----------
    mastorad : float
        Numerical factor to convert milliarcseconds ro radians.
    """

    def __init__(self):
        """
        Class constructor/initializer.
        """
        self.mastorad = np.pi / (180 * 3600 * 1000)

    def propagate_astrometry(
        self, phi, theta, parallax, muphistar, mutheta, vrad, t0, t1
    ):
        """
        Propagate the astrometric parameters of a source from the reference epoch t0 to
        the new epoch t1.

        Parameters
        ----------
        phi : float
            Longitude at reference epoch (radians).
        theta : float
            Latitude at reference epoch (radians).
        parallax : float
            Parallax at the reference epoch (mas).
        muphistar : float
            Proper motion in longitude (including np.cos(latitude) term) at reference
            epoch (mas/yr).
        mutheta : float
            Proper motion in latitude at reference epoch (mas/yr).
        vrad : float
            Radial velocity at reference epoch (km/s). Can be set to 0 km/s if not known.
        t0 : float
            Reference epoch (Julian years).
        t1 : float
            New epoch (Julian years).

        Returns
        -------
        phi1, theta1, parallax1, muphistar1, mutheta1, murad1 : float or array
            Astrometric parameters, including the "radial proper motion" (NOT the radial
            velocity), at the new epoch.
        """

        t = t1 - t0
        p0, q0, r0 = normal_triad(phi, theta)

        # Convert input data to units of radians and Julian year. Use ICRS coordinate
        # names internally to avoid errors in translating the formulae to code.
        pmra0 = muphistar * self.mastorad
        pmdec0 = mutheta * self.mastorad
        pmr0 = vrad * parallax / au_km_year_per_sec * self.mastorad
        pmtot0sqr = (muphistar**2 + mutheta**2) * self.mastorad**2

        # Proper motion vector
        pmvec0 = pmra0 * p0 + pmdec0 * q0

        f = (1 + 2 * pmr0 * t + (pmtot0sqr + pmr0**2) * t**2) ** (-0.5)
        u = (r0 * (1 + pmr0 * t) + pmvec0 * t) * f

        _, phi1, theta1 = cartesian_to_spherical(u[0], u[1], u[2])
        parallax1 = parallax * f
        pmr1 = (pmr0 + (pmtot0sqr + pmr0**2) * t) * f**2
        pmvec1 = (pmvec0 * (1 + pmr0 * t) - r0 * pmr0**2 * t) * f**3
        p1, q1, r1 = normal_triad(phi1, theta1)
        muphistar1 = np.sum(p1 * pmvec1 / self.mastorad, axis=0)
        mutheta1 = np.sum(q1 * pmvec1 / self.mastorad, axis=0)
        murad1 = pmr1 / self.mastorad

        return phi1, theta1, parallax1, muphistar1, mutheta1, murad1

    def propagate_pos(self, phi, theta, parallax, muphistar, mutheta, vrad, t0, t1):
        """
        Propagate the position of a source from the reference epoch t0 to the new epoch
        t1.

        Parameters
        ----------
        phi : float
            Longitude at reference epoch (radians).
        theta : float
            Latitude at reference epoch (radians).
        parallax : float
            Parallax at the reference epoch (mas).
        muphistar : float
            Proper motion in longitude (including np.cos(latitude) term) at reference
            epoch (mas/yr).
        mutheta : float
            Proper motion in latitude at reference epoch (mas/yr).
        vrad : float
            Radial velocity at reference epoch (km/s). Can be set to 0 km/s if not known.
        t0 : float
            Reference epoch (Julian years).
        t1 : float
            New epoch (Julian years).

        Returns
        -------
        phi, theta : float
            Coordinates the new epoch (in radians)
        """
        (phi1, theta1, _, _, _, _,) = self.propagate_astrometry(
            phi, theta, parallax, muphistar, mutheta, vrad, t0, t1
        )
        return phi1, theta1

    def propagate_astrometry_and_covariance_matrix(self, a0, c0, t0, t1):
        """
        Propagate the astrometric parameters iand radial proper motion, as well as the
        corresponding covariance matrix, from epoch t0 to epoch t1.

        Code based on the Hipparcos Fortran implementation by Lennart Lindegren.

        Parameters
        ----------
        a0 : array_like
            6-element vector: (phi, theta, parallax, muphistar, mutheta, vrad) in units
            of (radians, radians, mas, mas/yr, mas/yr, km/s). Shape of a should be (6,)
            or (6,N), with N the number of sources for which the astrometric parameters
            are provided. The value of vrad can be set to 0 km/s if the radial velocity is not know. An appropriate uncertainty should then be provided (see below).
        c0 : array_like
            Covariance matrix stored in a 6x6 element array. This can be constructed
            from the columns listed in the Gaia catalogue. The units are [mas^2,
            mas^2/yr, mas^2/yr^2] for the various elements. Note that the elements in
            the 6th row and column should be:
            c[6,i] = c[i,6] = c[i,3] * vrad / auKmYearPerSec
            for i = 1,..,5 and
            c[6,6] = c[3,3] * (vrad^2+vrad_error^2) / auKmYearPerSec^2 +
            (parallax*vrad_error/auKmYearPerSec)^2
            The shape of c0 should be (6,6) or (N,6,6). If the radial velocity is not know the uncertainty on the radial velocity (vrad_error) should be set to the velocity dispersion of the population the source is drawn from.
        t0 : float
            Reference epoch (Julian years).
        t1 : float
            New epoch (Julian years).

        Returns
        -------
        a, c : array
            Astrometric parameters, including the "radial proper motion" (NOT the radial
            velocity), and the covariance matrix at the new epoch as a 2D matrix with
            the new variances on the diagonal and the covariance in the off-diagonal
            elements.
        """

        zero, one, two, three = 0, 1, 2, 3
        tau = t1 - t0

        # Calculate the normal triad [p0 q0 r0] at t0
        p0, q0, r0 = normal_triad(a0[0], a0[1])

        # Convert to internal units (radians, Julian year)
        par0 = a0[2] * self.mastorad
        pma0 = a0[3] * self.mastorad
        pmd0 = a0[4] * self.mastorad
        pmr0 = a0[5] * a0[2] / au_km_year_per_sec * self.mastorad

        # Proper motion vector
        pmvec0 = pma0 * p0 + pmd0 * q0

        # Auxiliary quantities
        tau2 = tau * tau
        pm02 = pma0**2 + pmd0**2
        w = one + pmr0 * tau
        f2 = one / (one + two * pmr0 * tau + (pm02 + pmr0**2) * tau2)
        f = np.sqrt(f2)
        f3 = f2 * f
        f4 = f2 * f2

        # Position vector and parallax at t1
        u = (r0 * w + pmvec0 * tau) * f
        _, ra, dec = cartesian_to_spherical(u[0], u[1], u[2])
        par = par0 * f

        # Proper motion vector and radial proper motion at t1
        pmvec = (pmvec0 * (one + pmr0 * tau) - r0 * pmr0**2 * tau) * f3
        pmr = (pmr0 + (pm02 + pmr0**2) * tau) * f2

        # Normal triad at t1
        p, q, r = normal_triad(ra, dec)

        # Convert parameters at t1 to external units (mas, Julian year)
        pma = np.sum(p * pmvec, axis=0)
        pmd = np.sum(q * pmvec, axis=0)

        a = np.zeros_like(a0)
        a[0] = ra
        a[1] = dec
        a[2] = par / self.mastorad
        a[3] = pma / self.mastorad
        a[4] = pmd / self.mastorad
        a[5] = pmr / self.mastorad

        # Auxiliary quantities for the partial derivatives

        pmz = pmvec0 * f - three * pmvec * w
        pp0 = np.sum(p * p0, axis=0)
        pq0 = np.sum(p * q0, axis=0)
        pr0 = np.sum(p * r0, axis=0)
        qp0 = np.sum(q * p0, axis=0)
        qq0 = np.sum(q * q0, axis=0)
        qr0 = np.sum(q * r0, axis=0)
        ppmz = np.sum(p * pmz, axis=0)
        qpmz = np.sum(q * pmz, axis=0)

        jacobian = np.zeros_like(c0)
        if c0.ndim == 2:
            jacobian = jacobian[np.newaxis, :, :]

        # Partial derivatives
        jacobian[:, 0, 0] = pp0 * w * f - pr0 * pma0 * tau * f
        jacobian[:, 0, 1] = pq0 * w * f - pr0 * pmd0 * tau * f
        jacobian[:, 0, 2] = zero
        jacobian[:, 0, 3] = pp0 * tau * f
        jacobian[:, 0, 4] = pq0 * tau * f
        jacobian[:, 0, 5] = -pma * tau2

        jacobian[:, 1, 0] = qp0 * w * f - qr0 * pma0 * tau * f
        jacobian[:, 1, 1] = qq0 * w * f - qr0 * pmd0 * tau * f
        jacobian[:, 1, 2] = zero
        jacobian[:, 1, 3] = qp0 * tau * f
        jacobian[:, 1, 4] = qq0 * tau * f
        jacobian[:, 1, 5] = -pmd * tau2

        jacobian[:, 2, 0] = zero
        jacobian[:, 2, 1] = zero
        jacobian[:, 2, 2] = f
        jacobian[:, 2, 3] = -par * pma0 * tau2 * f2
        jacobian[:, 2, 4] = -par * pmd0 * tau2 * f2
        jacobian[:, 2, 5] = -par * w * tau * f2

        jacobian[:, 3, 0] = -pp0 * pm02 * tau * f3 - pr0 * pma0 * w * f3
        jacobian[:, 3, 1] = -pq0 * pm02 * tau * f3 - pr0 * pmd0 * w * f3
        jacobian[:, 3, 2] = zero
        jacobian[:, 3, 3] = (
            pp0 * w * f3 - two * pr0 * pma0 * tau * f3 - three * pma * pma0 * tau2 * f2
        )
        jacobian[:, 3, 4] = (
            pq0 * w * f3 - two * pr0 * pmd0 * tau * f3 - three * pma * pmd0 * tau2 * f2
        )
        jacobian[:, 3, 5] = ppmz * tau * f2

        jacobian[:, 4, 0] = -qp0 * pm02 * tau * f3 - qr0 * pma0 * w * f3
        jacobian[:, 4, 1] = -qq0 * pm02 * tau * f3 - qr0 * pmd0 * w * f3
        jacobian[:, 4, 2] = zero
        jacobian[:, 4, 3] = (
            qp0 * w * f3 - two * qr0 * pma0 * tau * f3 - three * pmd * pma0 * tau2 * f2
        )
        jacobian[:, 4, 4] = (
            qq0 * w * f3 - two * qr0 * pmd0 * tau * f3 - three * pmd * pmd0 * tau2 * f2
        )
        jacobian[:, 4, 5] = qpmz * tau * f2

        jacobian[:, 5, 0] = zero
        jacobian[:, 5, 1] = zero
        jacobian[:, 5, 2] = zero
        jacobian[:, 5, 3] = two * pma0 * w * tau * f4
        jacobian[:, 5, 4] = two * pmd0 * w * tau * f4
        jacobian[:, 5, 5] = (w**2 - pm02 * tau2) * f4

        jacobian_transposed = np.zeros_like(jacobian)
        for i in range(jacobian.shape[0]):
            jacobian_transposed[i] = jacobian[i].T

        if c0.ndim == 2:
            c = np.matmul(
                jacobian, np.matmul(c0[np.newaxis, :, :], jacobian_transposed)
            )
        else:
            c = np.matmul(jacobian, np.matmul(c0, jacobian_transposed))

        return a, np.squeeze(c)
