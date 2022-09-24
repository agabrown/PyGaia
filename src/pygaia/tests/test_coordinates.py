"""
Unit tests for the coordinates module.
"""
import numpy as np
from numpy.random import rand
import numpy.testing as npt

from pygaia.astrometry.constants import au_km_year_per_sec
from pygaia.astrometry.coordinates import CoordinateTransformation
from pygaia.astrometry.coordinates import (
    Transformations,
    EpochPropagation,
    angular_distance,
)
from pygaia.astrometry.vectorastrometry import (
    cartesian_to_spherical,
    astrometry_to_phase_space,
    phase_space_to_astrometry,
)
from pygaia.utils import construct_covariance_matrix

rng = np.random.default_rng()


class TestCoordinates(np.testing.TestCase):
    def setUp(self):
        self.basisVectorXValues = np.array([1, 0, 0])
        self.basisVectorYValues = np.array([0, 1, 0])
        self.basisVectorZValues = np.array([0, 0, 1])

        self.expectedMatIcrsToGal = np.array(
            [[-0.05, -0.87, -0.48], [+0.49, -0.44, +0.75], [-0.87, -0.20, +0.46]]
        )
        self.expectedMatEclipticToIcrs = np.array(
            [[1.0, 0.0, 0.0], [0.0, 0.92, -0.40], [0.0, 0.40, 0.92]]
        )
        self.expectedMatGalacticToEcliptic = np.array(
            [[-0.05, 0.49, -0.87], [-0.99, -0.11, 0.00], [-0.1, 0.86, 0.50]]
        )

    def test_galacticToIcrs(self):
        """
        Verify correctness of transformations from the Galactic to ICRS coordinate
        systems.
        """
        ct = CoordinateTransformation(Transformations.GAL2ICRS)
        x, y, z = ct.transform_cartesian_coordinates(
            self.basisVectorXValues, self.basisVectorYValues, self.basisVectorZValues
        )
        npt.assert_array_almost_equal(x, self.expectedMatIcrsToGal[:, 0], decimal=2)
        npt.assert_array_almost_equal(y, self.expectedMatIcrsToGal[:, 1], decimal=2)
        npt.assert_array_almost_equal(z, self.expectedMatIcrsToGal[:, 2], decimal=2)

        r, expectedAlpha, expectedDelta = cartesian_to_spherical(
            self.expectedMatIcrsToGal[:, 0],
            self.expectedMatIcrsToGal[:, 1],
            self.expectedMatIcrsToGal[:, 2],
        )
        alpha, delta = ct.transform_sky_coordinates(
            np.array([0.0, np.pi / 2.0, 0.0]), np.array([0.0, 0.0, np.pi / 2.0])
        )
        npt.assert_array_almost_equal(expectedAlpha, alpha, decimal=1)
        npt.assert_array_almost_equal(expectedDelta, delta, decimal=1)

        alpha, delta = ct.transform_sky_coordinates(0.0, np.pi / 2.0)
        npt.assert_allclose(alpha / np.pi * 180, 192.9, atol=1.0)
        npt.assert_allclose(delta / np.pi * 180, 27.1, atol=1.0)

    def test_icrsToGalactic(self):
        """
        Verify correctness of transformations from the ICRS to Galactic coordinate
        systems.
        """
        ct = CoordinateTransformation(Transformations.ICRS2GAL)
        x, y, z = ct.transform_cartesian_coordinates(
            self.basisVectorXValues, self.basisVectorYValues, self.basisVectorZValues
        )
        npt.assert_array_almost_equal(x, self.expectedMatIcrsToGal[0, :], decimal=2)
        npt.assert_array_almost_equal(y, self.expectedMatIcrsToGal[1, :], decimal=2)
        npt.assert_array_almost_equal(z, self.expectedMatIcrsToGal[2, :], decimal=2)

        r, expectedGalon, expectedGalat = cartesian_to_spherical(
            self.expectedMatIcrsToGal[0, :],
            self.expectedMatIcrsToGal[1, :],
            self.expectedMatIcrsToGal[2, :],
        )
        galon, galat = ct.transform_sky_coordinates(
            np.array([0.0, np.pi / 2.0, 0.0]), np.array([0.0, 0.0, np.pi / 2.0])
        )
        npt.assert_array_almost_equal(expectedGalon, galon, decimal=1)
        npt.assert_array_almost_equal(expectedGalat, galat, decimal=1)

    def test_icrsToEcliptic(self):
        """
        Verify correctness of transformations from the ICRS to Ecliptic coordinate
        systems.
        """
        ct = CoordinateTransformation(Transformations.ICRS2ECL)
        x, y, z = ct.transform_cartesian_coordinates(
            self.basisVectorXValues, self.basisVectorYValues, self.basisVectorZValues
        )
        npt.assert_array_almost_equal(
            x, self.expectedMatEclipticToIcrs[:, 0], decimal=2
        )
        npt.assert_array_almost_equal(
            y, self.expectedMatEclipticToIcrs[:, 1], decimal=2
        )
        npt.assert_array_almost_equal(
            z, self.expectedMatEclipticToIcrs[:, 2], decimal=2
        )

        r, expectedLambda, expectedBeta = cartesian_to_spherical(
            self.expectedMatEclipticToIcrs[:, 0],
            self.expectedMatEclipticToIcrs[:, 1],
            self.expectedMatEclipticToIcrs[:, 2],
        )
        lambdaEcl, betaEcl = ct.transform_sky_coordinates(
            np.array([0.0, np.pi / 2.0, 0.0]), np.array([0.0, 0.0, np.pi / 2.0])
        )
        npt.assert_array_almost_equal(expectedLambda, lambdaEcl, decimal=1)
        npt.assert_array_almost_equal(expectedBeta, betaEcl, decimal=1)

        alpha = 2.0 * np.pi * rng.uniform(size=100)
        delta = -np.pi / 2.0 + np.pi * rng.uniform(size=100)
        lambdaEcl, betaEcl = ct.transform_sky_coordinates(alpha, delta)
        result = 0.9175 * np.sin(delta) - 0.3978 * np.sin(alpha) * np.cos(delta)
        npt.assert_array_almost_equal(result, np.sin(betaEcl), decimal=2)

    def test_eclipticToIcrs(self):
        """
        Verify correctness of transformations from the Ecliptic to ICRS coordinate
        systems.
        """
        ct = CoordinateTransformation(Transformations.ECL2ICRS)
        x, y, z = ct.transform_cartesian_coordinates(
            self.basisVectorXValues, self.basisVectorYValues, self.basisVectorZValues
        )
        npt.assert_array_almost_equal(
            x, self.expectedMatEclipticToIcrs[0, :], decimal=2
        )
        npt.assert_array_almost_equal(
            y, self.expectedMatEclipticToIcrs[1, :], decimal=2
        )
        npt.assert_array_almost_equal(
            z, self.expectedMatEclipticToIcrs[2, :], decimal=2
        )

        r, expectedAlpha, expectedDelta = cartesian_to_spherical(
            self.expectedMatEclipticToIcrs[0, :],
            self.expectedMatEclipticToIcrs[1, :],
            self.expectedMatEclipticToIcrs[2, :],
        )
        alpha, delta = ct.transform_sky_coordinates(
            np.array([0.0, np.pi / 2.0, 0.0]), np.array([0.0, 0.0, np.pi / 2.0])
        )
        npt.assert_array_almost_equal(expectedAlpha, alpha, decimal=1)
        npt.assert_array_almost_equal(expectedDelta, delta, decimal=1)

        lambdaEcl = 2.0 * np.pi * rng.uniform(size=100)
        betaEcl = -np.pi / 2.0 + np.pi * rng.uniform(size=100)
        alpha, delta = ct.transform_sky_coordinates(lambdaEcl, betaEcl)
        result = 0.9175 * np.sin(betaEcl) + 0.3978 * np.sin(lambdaEcl) * np.cos(betaEcl)
        npt.assert_array_almost_equal(result, np.sin(delta), decimal=2)

    def test_eclipticToGalactic(self):
        """
        Verify correctness of transformations from the Ecliptic to Galactic coordinate
        systems.
        """
        ct = CoordinateTransformation(Transformations.ECL2GAL)
        x, y, z = ct.transform_cartesian_coordinates(
            self.basisVectorXValues, self.basisVectorYValues, self.basisVectorZValues
        )
        npt.assert_array_almost_equal(
            x, self.expectedMatGalacticToEcliptic[:, 0], decimal=2
        )
        npt.assert_array_almost_equal(
            y, self.expectedMatGalacticToEcliptic[:, 1], decimal=2
        )
        npt.assert_array_almost_equal(
            z, self.expectedMatGalacticToEcliptic[:, 2], decimal=2
        )

        r, expectedGalon, expectedGalat = cartesian_to_spherical(
            self.expectedMatGalacticToEcliptic[:, 0],
            self.expectedMatGalacticToEcliptic[:, 1],
            self.expectedMatGalacticToEcliptic[:, 2],
        )
        galon, galat = ct.transform_sky_coordinates(
            np.array([0.0, np.pi / 2.0, 0.0]), np.array([0.0, 0.0, np.pi / 2.0])
        )
        npt.assert_array_almost_equal(expectedGalon, galon, decimal=1)
        npt.assert_array_almost_equal(expectedGalat, galat, decimal=1)

    def test_galacticToEcliptic(self):
        """
        Verify correctness of transformations from the Galactic to Ecliptic coordinate
        systems.
        """
        ct = CoordinateTransformation(Transformations.GAL2ECL)
        x, y, z = ct.transform_cartesian_coordinates(
            self.basisVectorXValues, self.basisVectorYValues, self.basisVectorZValues
        )
        npt.assert_array_almost_equal(
            x, self.expectedMatGalacticToEcliptic[0, :], decimal=2
        )
        npt.assert_array_almost_equal(
            y, self.expectedMatGalacticToEcliptic[1, :], decimal=2
        )
        npt.assert_array_almost_equal(
            z, self.expectedMatGalacticToEcliptic[2, :], decimal=2
        )

        r, expectedLambda, expectedBeta = cartesian_to_spherical(
            self.expectedMatGalacticToEcliptic[0, :],
            self.expectedMatGalacticToEcliptic[1, :],
            self.expectedMatGalacticToEcliptic[2, :],
        )
        lambdaEcl, betaEcl = ct.transform_sky_coordinates(
            np.array([0.0, np.pi / 2.0, 0.0]), np.array([0.0, 0.0, np.pi / 2.0])
        )
        #
        # Note the abs() in the line below is required to avoid an error when -pi and pi
        # are compared. The better solution of course is to write an npt.assert
        # function that can handle modulo 2*pi cases.
        #
        npt.assert_array_almost_equal(abs(expectedLambda), abs(lambdaEcl), decimal=1)
        npt.assert_array_almost_equal(expectedBeta, betaEcl, decimal=1)

        galon = 2.0 * np.pi * rng.uniform(size=100)
        galat = -np.pi / 2.0 + np.pi * rng.uniform(size=100)
        lambdaEcl, betaEcl = ct.transform_sky_coordinates(galon, galat)
        phi = 6.38 / 180.0 * np.pi
        result = 0.4971 * np.sin(galat) + 0.8677 * np.sin(galon - phi) * np.cos(galat)
        npt.assert_array_almost_equal(result, np.sin(betaEcl), decimal=2)

    def test_transformProperMotions(self):
        """
        Verify the correct implementation of the direct transformation of proper
        motions.
        """
        ct = CoordinateTransformation(Transformations.ICRS2GAL)
        nTests = 100
        phi = 2.0 * np.pi * rng.uniform(size=nTests)
        theta = -np.pi / 2.0 + np.pi * rng.uniform(size=nTests)
        parallax = rng.uniform(size=nTests) * 99.0 + 1.0
        muphistar = -30.0 + 60.0 * rng.uniform(size=nTests)
        mutheta = -30.0 + 60.0 * rng.uniform(size=nTests)
        vrad = -200.0 + 400.0 * rng.uniform(size=nTests)

        x, y, z, vx, vy, vz = astrometry_to_phase_space(
            phi, theta, parallax, muphistar, mutheta, vrad
        )
        xrot, yrot, zrot = ct.transform_cartesian_coordinates(x, y, z)
        vxrot, vyrot, vzrot = ct.transform_cartesian_coordinates(vx, vy, vz)
        (
            phiRot,
            thetaRot,
            parRot,
            muphistarRotExpected,
            muthetaRotExpected,
            vradRot,
        ) = phase_space_to_astrometry(xrot, yrot, zrot, vxrot, vyrot, vzrot)

        muphistarRot, muthetaRot = ct.transform_proper_motions(
            phi, theta, muphistar, mutheta
        )

        npt.assert_array_almost_equal(muphistarRotExpected, muphistarRot, decimal=2)
        npt.assert_array_almost_equal(muthetaRotExpected, muthetaRot, decimal=2)

        #
        # Test scalar version of method.
        #
        for i in range(nTests):
            muphistarRot, muthetaRot = ct.transform_proper_motions(
                phi[i], theta[i], muphistar[i], mutheta[i]
            )
            npt.assert_almost_equal(muphistarRotExpected[i], muphistarRot, decimal=2)
            npt.assert_almost_equal(muthetaRotExpected[i], muthetaRot, decimal=2)

    def test_transformSkyCoordinateErrors(self):
        """
        Verify that the transformed covariance matrix for the positions remains a
        covariance matrix.
        """
        ct = CoordinateTransformation(Transformations.ICRS2GAL)
        nTests = 100
        phi = 2.0 * np.pi * rng.uniform(size=nTests)
        theta = -np.pi / 2.0 + np.pi * rng.uniform(size=nTests)
        sigPhiStar, sigTheta, rhoPhiTheta = self._generateRandomCovarianceMatrices(
            nTests
        )

        sigPhiStarRot, sigThetaRot, rhoPhiThetaRot = ct.transform_sky_coordinate_errors(
            phi, theta, sigPhiStar, sigTheta, rho_phi_theta=rhoPhiTheta
        )
        for i in range(nTests):
            self.assertGreater(sigPhiStarRot[i], 0.0)
            self.assertGreater(sigThetaRot[i], 0.0)
            self.assertTrue(-1 <= rhoPhiThetaRot[i] <= 1)

        sigPhiStarRot, sigThetaRot, rhoPhiThetaRot = ct.transform_sky_coordinate_errors(
            phi, theta, sigPhiStar, sigTheta
        )
        for i in range(nTests):
            self.assertGreater(sigPhiStarRot[i], 0.0)
            self.assertGreater(sigThetaRot[i], 0.0)
            self.assertTrue(-1 <= rhoPhiThetaRot[i] <= 1)

        sigPhiStarRot, sigThetaRot, rhoPhiThetaRot = ct.transform_sky_coordinate_errors(
            phi, theta, sigPhiStar, sigTheta, rho_phi_theta=0.7
        )
        for i in range(nTests):
            self.assertGreater(sigPhiStarRot[i], 0.0)
            self.assertGreater(sigThetaRot[i], 0.0)
            self.assertTrue(-1 <= rhoPhiThetaRot[i] <= 1)

        for i in range(nTests):
            (
                sigPhiStarRot,
                sigThetaRot,
                rhoPhiThetaRot,
            ) = ct.transform_sky_coordinate_errors(
                phi[i], theta[i], sigPhiStar[i], sigTheta[i], rho_phi_theta=0.7
            )
            self.assertGreater(sigPhiStarRot, 0.0)
            self.assertGreater(sigThetaRot, 0.0)
            self.assertTrue(-1 <= rhoPhiThetaRot <= 1)

    def test_transformProperMotionErrors(self):
        """
        Verify that the transformed covariance matrix for the proper motions remains a
        covariance matrix.
        """
        ct = CoordinateTransformation(Transformations.ICRS2GAL)
        nTests = 100
        phi = 2.0 * np.pi * rng.uniform(size=nTests)
        theta = -np.pi / 2.0 + np.pi * rng.uniform(size=nTests)
        sigPhiStar, sigTheta, rhoPhiTheta = self._generateRandomCovarianceMatrices(
            nTests
        )

        sigPhiStarRot, sigThetaRot, rhoPhiThetaRot = ct.transform_proper_motion_errors(
            phi, theta, sigPhiStar, sigTheta, rho_muphi_mutheta=rhoPhiTheta
        )
        for i in range(nTests):
            self.assertGreater(sigPhiStarRot[i], 0.0)
            self.assertGreater(sigThetaRot[i], 0.0)
            self.assertTrue(-1 <= rhoPhiTheta[i] <= 1)

        sigPhiStarRot, sigThetaRot, rhoPhiThetaRot = ct.transform_proper_motion_errors(
            phi, theta, sigPhiStar, sigTheta
        )
        for i in range(nTests):
            self.assertGreater(sigPhiStarRot[i], 0.0)
            self.assertGreater(sigThetaRot[i], 0.0)
            self.assertTrue(-1 <= rhoPhiTheta[i] <= 1)

    def test_propagate_astrometry(self):
        """
        Verify correct working of the epoch propagation class.
        """
        radtomas = 180 * 3600 * 1000 / np.pi

        ep = EpochPropagation()
        phi = 0.0
        theta = 0.0
        parallax = 10.0
        muphistar = 1.0
        mutheta = 0.0
        vrad = 0.0
        t0 = 2015.5
        t1 = 2020.5
        t2 = 2010.0
        phi1, theta1, parallax1, muphistar1, mutheta1, pmr1 = ep.propagate_astrometry(
            phi, theta, parallax, muphistar, mutheta, vrad, t0, t1
        )
        self.assertGreater(phi1, 0.0)
        self.assertEquals(theta1, 0.0)
        self.assertLess(parallax1, parallax)
        self.assertLess(muphistar1, muphistar)
        self.assertEquals(mutheta1, 0.0)
        self.assertGreater(pmr1, 0.0)

        # For propagation into the future check that in absence of radial motion the
        # parallaxes and proper motions always diminish with respect to the reference
        # epoch, while the radial proper motion always becomes positive non-zero.
        nTests = 100
        phi = 2.0 * np.pi * rng.uniform(size=nTests)
        theta = np.arcsin(-1 + 2 * rng.uniform(size=nTests))
        parallax = rng.uniform(size=nTests) * 20.0 + 1.0
        muphistar = rng.uniform(size=nTests) * 5.0 - 10.0
        mutheta = rng.uniform(size=nTests) * 5.0 - 10.0
        vrad = np.zeros_like(phi)
        phi1, theta1, parallax1, muphistar1, mutheta1, pmr1 = ep.propagate_astrometry(
            phi, theta, parallax, muphistar, mutheta, vrad, t0, t1
        )
        npt.assert_array_less(parallax1, parallax)
        npt.assert_array_less(
            muphistar1**2 + mutheta1**2, muphistar**2 + mutheta**2
        )
        npt.assert_array_less(0.0, pmr1)

        # For propagation into the past check that in absence of radial motion the
        # parallaxes and proper motions always diminish with respect to the reference
        # epoch, while the radial proper motion always becomes negative non-zero.
        phi2, theta2, parallax2, muphistar2, mutheta2, pmr2 = ep.propagate_astrometry(
            phi, theta, parallax, muphistar, mutheta, vrad, t0, t2
        )
        npt.assert_array_less(parallax2, parallax)
        npt.assert_array_less(
            muphistar2**2 + mutheta2**2, muphistar**2 + mutheta**2
        )
        npt.assert_array_less(pmr2, 0.0)

        # For propagation into the future check that In the absence of proper motion and
        # for positive radial velocity and positive parallax at the current epoch: the
        # parallax always decreases, the radial proper motion decreases.
        muphistar = np.zeros_like(phi)
        mutheta = np.zeros_like(phi)
        vrad = rng.uniform(size=nTests) * 50.0 + 1.0
        pmr = vrad * parallax / au_km_year_per_sec
        phi1, theta1, parallax1, muphistar1, mutheta1, pmr1 = ep.propagate_astrometry(
            phi, theta, parallax, muphistar, mutheta, vrad, t0, t1
        )
        npt.assert_almost_equal(phi1, phi, 12)
        npt.assert_almost_equal(theta1, theta, 12)
        npt.assert_array_less(parallax1, parallax)
        npt.assert_almost_equal(muphistar1, muphistar, 12)
        npt.assert_almost_equal(mutheta1, mutheta, 12)
        npt.assert_array_less(pmr1, pmr)

        # For propagation into the past check that In the absence of proper motion and
        # for positive radial velocity and positive parallax at the current epoch: the
        # parallax always increases, the radial proper motion always increases.
        phi2, theta2, parallax2, muphistar2, mutheta2, pmr2 = ep.propagate_astrometry(
            phi, theta, parallax, muphistar, mutheta, vrad, t0, t2
        )
        npt.assert_almost_equal(phi2, phi, 12)
        npt.assert_almost_equal(theta2, theta, 12)
        npt.assert_array_less(parallax, parallax2)
        npt.assert_almost_equal(muphistar2, muphistar, 12)
        npt.assert_almost_equal(mutheta2, mutheta, 12)
        npt.assert_array_less(pmr, pmr2)

        # For propagation into the future check that In the absence of proper motion and
        # for negative radial velocity and positive parallax at the current epoch: the
        # parallax always increases, the radial proper motion always decreases.
        vrad = -1.0 - 70 * rng.uniform(size=nTests)
        pmr = vrad * parallax / au_km_year_per_sec
        phi1, theta1, parallax1, muphistar1, mutheta1, pmr1 = ep.propagate_astrometry(
            phi, theta, parallax, muphistar, mutheta, vrad, t0, t1
        )
        npt.assert_almost_equal(phi1, phi, 12)
        npt.assert_almost_equal(theta1, theta, 12)
        npt.assert_array_less(parallax, parallax1)
        npt.assert_almost_equal(muphistar1, muphistar, 12)
        npt.assert_almost_equal(mutheta1, mutheta, 12)
        npt.assert_array_less(pmr1, pmr)

        # For propagation into the past check that In the absence of proper motion and
        # for negative radial velocity and positive parallax at the current epoch: the
        # parallax always decreases, and the radial proper motion always increases.
        phi2, theta2, parallax2, muphistar2, mutheta2, pmr2 = ep.propagate_astrometry(
            phi, theta, parallax, muphistar, mutheta, vrad, t0, t2
        )
        npt.assert_almost_equal(phi2, phi, 12)
        npt.assert_almost_equal(theta2, theta, 12)
        npt.assert_array_less(parallax2, parallax)
        npt.assert_almost_equal(muphistar2, muphistar, 12)
        npt.assert_almost_equal(mutheta2, mutheta, 12)
        npt.assert_array_less(pmr, pmr2)

        # Check that angular distance between old and new position is close to that
        # obtained from a naive epoch propagation.
        nTests = 10000
        phi = 2.0 * np.pi * rng.uniform(size=nTests)
        theta = np.arcsin(-1 + 2 * rng.uniform(size=nTests))
        parallax = rng.uniform(size=nTests) * 20.0
        muphistar = rng.uniform(size=nTests) * 200 - 100
        mutheta = rng.uniform(size=nTests) * 200 - 100
        vrad = rng.uniform(size=nTests) * 400 - 200
        phi1, theta1, parallax1, muphistar1, mutheta1, pmr1 = ep.propagate_astrometry(
            phi, theta, parallax, muphistar, mutheta, vrad, t0, t1
        )
        phi2, theta2, parallax2, muphistar2, mutheta2, pmr2 = ep.propagate_astrometry(
            phi, theta, parallax, muphistar, mutheta, vrad, t0, t2
        )
        rho1 = angular_distance(phi, theta, phi1, theta1) * radtomas
        rho2 = angular_distance(phi, theta, phi2, theta2) * radtomas
        rho_naive1 = np.sqrt(muphistar**2 + mutheta**2) * abs(t1 - t0)
        rho_naive2 = np.sqrt(muphistar**2 + mutheta**2) * abs(t2 - t0)
        npt.assert_allclose(rho_naive1, rho1, rtol=0.1)
        npt.assert_allclose(rho_naive2, rho2, rtol=0.1)

        # Check that for zero parallax the propgation still works, leaving parallax
        # unchanged.
        nTests = 100
        phi = 2.0 * np.pi * rng.uniform(size=nTests)
        theta = np.arcsin(-1 + 2 * rng.uniform(size=nTests))
        parallax = np.zeros_like(phi)
        muphistar = rng.uniform(size=nTests) * 2 - 1
        mutheta = rng.uniform(size=nTests) * 2 - 1
        vrad = rng.uniform(size=nTests) * 40 - 20
        phi1, theta1, parallax1, muphistar1, mutheta1, pmr1 = ep.propagate_astrometry(
            phi, theta, parallax, muphistar, mutheta, vrad, t0, t1
        )
        phi2, theta2, parallax2, muphistar2, mutheta2, pmr2 = ep.propagate_astrometry(
            phi, theta, parallax, muphistar, mutheta, vrad, t0, t2
        )
        npt.assert_equal(parallax1, 0.0)
        npt.assert_equal(parallax2, 0.0)
        npt.assert_array_less(0.0, pmr1)
        npt.assert_array_less(pmr2, 0.0)

        # For propagation into the future check that In the absence of proper motion and
        # for positive radial velocity and negative parallax at the current epoch: the
        # parallax always decreases, the radial proper motion decreases.
        nTests = 100
        phi = 2.0 * np.pi * rng.uniform(size=nTests)
        theta = np.arcsin(-1 + 2 * rng.uniform(size=nTests))
        parallax = -rng.uniform(size=nTests) * 20.0 - 1.0
        muphistar = np.zeros_like(phi)
        mutheta = np.zeros_like(phi)
        vrad = rng.uniform(size=nTests) * 50.0 + 1.0
        pmr = vrad * parallax / au_km_year_per_sec
        phi1, theta1, parallax1, muphistar1, mutheta1, pmr1 = ep.propagate_astrometry(
            phi, theta, parallax, muphistar, mutheta, vrad, t0, t1
        )
        npt.assert_almost_equal(phi1, phi, 12)
        npt.assert_almost_equal(theta1, theta, 12)
        npt.assert_array_less(parallax1, parallax)
        npt.assert_almost_equal(muphistar1, muphistar, 12)
        npt.assert_almost_equal(mutheta1, mutheta, 12)
        npt.assert_array_less(pmr1, pmr)

        # For propagation into the past check that In the absence of proper motion and
        # for positive radial velocity and negative parallax at the current epoch: the
        # parallax always increases, the radial proper motion always increases.
        phi2, theta2, parallax2, muphistar2, mutheta2, pmr2 = ep.propagate_astrometry(
            phi, theta, parallax, muphistar, mutheta, vrad, t0, t2
        )
        npt.assert_almost_equal(phi2, phi, 12)
        npt.assert_almost_equal(theta2, theta, 12)
        npt.assert_array_less(parallax, parallax2)
        npt.assert_almost_equal(muphistar2, muphistar, 12)
        npt.assert_almost_equal(mutheta2, mutheta, 12)
        npt.assert_array_less(pmr, pmr2)

        # For propagation into the future check that In the absence of proper motion and
        # for negative radial velocity and negative parallax at the current epoch: the
        # parallax always increases, the radial proper motion always decreases.
        vrad = -1.0 - 70 * rng.uniform(size=nTests)
        pmr = vrad * parallax / au_km_year_per_sec
        phi1, theta1, parallax1, muphistar1, mutheta1, pmr1 = ep.propagate_astrometry(
            phi, theta, parallax, muphistar, mutheta, vrad, t0, t1
        )
        npt.assert_almost_equal(phi1, phi, 12)
        npt.assert_almost_equal(theta1, theta, 12)
        npt.assert_array_less(parallax, parallax1)
        npt.assert_almost_equal(muphistar1, muphistar, 12)
        npt.assert_almost_equal(mutheta1, mutheta, 12)
        npt.assert_array_less(pmr1, pmr)

        # For propagation into the past check that In the absence of proper motion and
        # for negative radial velocity and negative parallax at the current epoch: the
        # parallax always decreases, and the radial proper motion always increases.
        phi2, theta2, parallax2, muphistar2, mutheta2, pmr2 = ep.propagate_astrometry(
            phi, theta, parallax, muphistar, mutheta, vrad, t0, t2
        )
        npt.assert_almost_equal(phi2, phi, 12)
        npt.assert_almost_equal(theta2, theta, 12)
        npt.assert_array_less(parallax2, parallax)
        npt.assert_almost_equal(muphistar2, muphistar, 12)
        npt.assert_almost_equal(mutheta2, mutheta, 12)
        npt.assert_array_less(pmr, pmr2)

    def test_propagate_astrometry_and_covariance_matrix(self):
        """
        Verify correct working of propagate_astrometry_and_covariance_matrix().
        """

        ep = EpochPropagation()
        t0 = 2015.5
        t1 = 2020.5
        t2 = 2010.0

        nTests = 10000
        a0 = np.zeros((6, nTests))
        a0[0] = 2.0 * np.pi * rng.uniform(size=nTests)
        a0[1] = np.arcsin(-1 + 2 * rng.uniform(size=nTests))
        a0[2] = rng.uniform(size=nTests) * 20.0
        a0[3] = rng.uniform(size=nTests) * 200 - 100
        a0[4] = rng.uniform(size=nTests) * 200 - 100
        a0[5] = rng.uniform(size=nTests) * 400 - 200

        covmat = np.array(
            [1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        )
        vrad_error = 2.0
        c0 = construct_covariance_matrix(
            np.tile(covmat, (nTests, covmat.size)),
            a0[2],
            a0[5],
            np.tile(vrad_error, nTests),
        )

        for i in range(nTests):
            self.assertTrue(self._is_pos_def(c0[i]))

        a1, c1 = ep.propagate_astrometry_and_covariance_matrix(a0, c0, t0, t1)
        a2, c2 = ep.propagate_astrometry_and_covariance_matrix(a0, c0, t0, t2)
        phi1, theta1, parallax1, muphistar1, mutheta1, pmr1 = ep.propagate_astrometry(
            a0[0], a0[1], a0[2], a0[3], a0[4], a0[5], t0, t1
        )
        phi2, theta2, parallax2, muphistar2, mutheta2, pmr2 = ep.propagate_astrometry(
            a0[0], a0[1], a0[2], a0[3], a0[4], a0[5], t0, t2
        )

        npt.assert_almost_equal(phi1, a1[0], 12)
        npt.assert_almost_equal(theta1, a1[1], 12)
        npt.assert_almost_equal(parallax1, a1[2], 12)
        npt.assert_almost_equal(muphistar1, a1[3], 12)
        npt.assert_almost_equal(mutheta1, a1[4], 12)
        npt.assert_almost_equal(pmr1, a1[5], 12)

        npt.assert_almost_equal(phi2, a2[0], 12)
        npt.assert_almost_equal(theta2, a2[1], 12)
        npt.assert_almost_equal(parallax2, a2[2], 12)
        npt.assert_almost_equal(muphistar2, a2[3], 12)
        npt.assert_almost_equal(mutheta2, a2[4], 12)
        npt.assert_almost_equal(pmr2, a2[5], 12)

        for i in range(nTests):
            self.assertTrue(self._is_pos_def(c1[i]))
            self.assertTrue(self._is_pos_def(c2[i]))

        a1, c1 = ep.propagate_astrometry_and_covariance_matrix(a0[:, 0], c0[0], t0, t1)
        a2, c2 = ep.propagate_astrometry_and_covariance_matrix(a0[:, 0], c0[0], t0, t2)
        phi1, theta1, parallax1, muphistar1, mutheta1, pmr1 = ep.propagate_astrometry(
            a0[0, 0], a0[1, 0], a0[2, 0], a0[3, 0], a0[4, 0], a0[5, 0], t0, t1
        )
        phi2, theta2, parallax2, muphistar2, mutheta2, pmr2 = ep.propagate_astrometry(
            a0[0, 0], a0[1, 0], a0[2, 0], a0[3, 0], a0[4, 0], a0[5, 0], t0, t2
        )

        npt.assert_almost_equal(phi1, a1[0], 12)
        npt.assert_almost_equal(theta1, a1[1], 12)
        npt.assert_almost_equal(parallax1, a1[2], 12)
        npt.assert_almost_equal(muphistar1, a1[3], 12)
        npt.assert_almost_equal(mutheta1, a1[4], 12)
        npt.assert_almost_equal(pmr1, a1[5], 12)

        npt.assert_almost_equal(phi2, a2[0], 12)
        npt.assert_almost_equal(theta2, a2[1], 12)
        npt.assert_almost_equal(parallax2, a2[2], 12)
        npt.assert_almost_equal(muphistar2, a2[3], 12)
        npt.assert_almost_equal(mutheta2, a2[4], 12)
        npt.assert_almost_equal(pmr2, a2[5], 12)

        self.assertTrue(self._is_pos_def(c1))
        self.assertTrue(self._is_pos_def(c2))

    def test_angular_distance(self):
        """
        Verify correct working of angular_distance function.

        Returns
        -------
        Nothing
        """
        ra1 = np.deg2rad([0, 0, 0, 35, 35, 65, 200, 80])
        dec1 = np.deg2rad([0, 0, 0, 0, 0, 30, 90, 0])
        ra2 = np.deg2rad([0, 0, 90, 35, 215, 245, 200, 350])
        dec2 = np.deg2rad([0, 90, 0, 10, 0, -30, 80, 0])

        rho, theta = angular_distance(ra1, dec1, ra2, dec2, return_posangle=True)
        expected_rho = np.deg2rad([0, 90, 90, 10, 180, 180, 10, 90])
        expected_theta = np.deg2rad([0, 0, 90, 0, 0, 0, 180, 270])

        npt.assert_array_almost_equal(rho, expected_rho, 12)
        npt.assert_array_almost_equal(theta, expected_theta, 12)

    @staticmethod
    def _generateRandomCovarianceMatrices(num):
        """
        Generate random covariance matrices through the transformation of diagonal
        positive definite matrices.

        Parameters
        ----------
        num : int
            Number of matrices to generate.

        Returns
        -------
        sigma1 : float
            Square root of variance of first variable.
        sigma2 : float
            Square root of variance of second variable.
        rho12 : float
            Correlation coefficient between errors on variables 1 and 2
        """
        varDiag1 = rng.uniform(size=num) * 99.0 + 1.0
        varDiag2 = rng.uniform(size=num) * 99.0 + 1.0
        rotAngle = 2.0 * np.pi * rng.uniform(size=num)
        sigma1 = np.zeros(num)
        sigma2 = np.zeros(num)
        rho12 = np.zeros(num)
        for i in range(num):
            rotationMatrix = np.array(
                [
                    [np.cos(rotAngle[i]), np.sin(rotAngle[i])],
                    [-np.sin(rotAngle[i]), np.cos(rotAngle[i])],
                ]
            )
            covMatrix = np.dot(
                rotationMatrix,
                np.dot(
                    np.diag([varDiag1[i], varDiag2[i]]), np.transpose(rotationMatrix)
                ),
            )
            sigma1[i] = np.sqrt(covMatrix[0, 0])
            sigma2[i] = np.sqrt(covMatrix[1, 1])
            rho12[i] = covMatrix[0, 1] / (sigma1[i] * sigma2[i])

        return sigma1, sigma2, rho12

    @staticmethod
    def _is_pos_def(A):
        """
        Check if matrix A is positive definite. Code from
        https://stackoverflow.com/questions/16266720/find-out-if-matrix-is-positive-definite-with-numpy
        """
        if np.allclose(A, A.T, rtol=1e-10, atol=1e-12):
            try:
                np.linalg.cholesky(A)
                return True
            except np.linalg.LinAlgError:
                return False
        else:
            return False
