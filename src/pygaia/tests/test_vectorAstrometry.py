"""
Unit tests for the vectorastrometry module.
"""

import numpy as np

from pygaia.astrometry.vectorastrometry import cartesian_to_spherical
from pygaia.astrometry.vectorastrometry import spherical_to_cartesian
from pygaia.astrometry.vectorastrometry import elementary_rotation_matrix
from pygaia.astrometry.vectorastrometry import normal_triad
from pygaia.astrometry.vectorastrometry import phase_space_to_astrometry
from pygaia.astrometry.vectorastrometry import astrometry_to_phase_space

rng = np.random.default_rng()


class test_vectorAstrometry(np.testing.TestCase):
    def test_cartesianToSpherical(self):
        """
        Check correct working of the Cartesian to Spherical coordinates transformation.
        """
        r, phi, theta = cartesian_to_spherical(
            np.array([1, 0, 0]), np.array([0, 1, 0]), np.array([0, 0, 1])
        )
        np.testing.assert_array_almost_equal(r, np.array([1.0, 1.0, 1.0]))
        np.testing.assert_array_almost_equal(phi, np.array([0.0, np.pi / 2.0, 0.0]))
        np.testing.assert_array_almost_equal(theta, np.array([0.0, 0.0, np.pi / 2.0]))

        r, phi, theta = cartesian_to_spherical(1.0, 1.0, 0.0)
        np.testing.assert_almost_equal(r, np.sqrt(2.0))
        np.testing.assert_almost_equal(phi, np.pi / 4.0)
        np.testing.assert_almost_equal(theta, 0.0)

        r, phi, theta = cartesian_to_spherical(1.0, 1.0, 1.0)
        np.testing.assert_almost_equal(r, np.sqrt(3.0))
        np.testing.assert_almost_equal(phi, np.pi / 4.0)
        np.testing.assert_almost_equal(theta, np.arcsin(1.0 / np.sqrt(3.0)))

        r, phi, theta = cartesian_to_spherical(1.0, 1.0, -1.0)
        np.testing.assert_almost_equal(r, np.sqrt(3.0))
        np.testing.assert_almost_equal(phi, np.pi / 4.0)
        np.testing.assert_almost_equal(theta, -np.arcsin(1.0 / np.sqrt(3.0)))

        r, phi, theta = cartesian_to_spherical(-1.0, -1.0, -1.0)
        np.testing.assert_almost_equal(r, np.sqrt(3.0))
        np.testing.assert_almost_equal(phi, np.pi / 4.0 + np.pi)
        np.testing.assert_almost_equal(theta, -np.arcsin(1.0 / np.sqrt(3.0)))

        np.testing.assert_raises(Exception, cartesian_to_spherical, 0.0, 0.0, 0.0)
        np.testing.assert_raises(
            Exception,
            cartesian_to_spherical,
            np.array([1, 0, 0, 0]),
            np.array([0, 1, 0, 0]),
            np.array([0, 0, 0, 1]),
        )

    def test_sphericalToCartesian(self):
        """
        Check correct working of the Spherical to Cartesian coordinates transformation.
        """
        x, y, z = spherical_to_cartesian(
            np.array([1.0, 2.0, 3.0, 5.0]),
            np.array([0.0, np.pi / 2.0, 0.0, 3.0 * np.pi / 4.0]),
            np.array([0.0, 0.0, np.pi / 2.0, np.pi / 6.0]),
        )
        np.testing.assert_array_almost_equal(
            x, np.array([1.0, 0.0, 0.0, -5.0 * np.sqrt(6.0) / 4.0])
        )
        np.testing.assert_array_almost_equal(
            y, np.array([0.0, 2.0, 0.0, 5.0 * np.sqrt(6.0) / 4.0])
        )
        np.testing.assert_array_almost_equal(z, np.array([0.0, 0.0, 3.0, 2.5]))

        x, y, z = spherical_to_cartesian(3.0, -np.pi / 6.0, -np.pi / 4.0)
        np.testing.assert_almost_equal(x, 3.0 * np.sqrt(6.0) / 4.0)
        np.testing.assert_almost_equal(y, -3.0 * np.sqrt(2.0) / 4.0)
        np.testing.assert_almost_equal(z, -3.0 * np.sqrt(2.0) / 2.0)

    def test_elementaryRotationMatrix(self):
        """
        Check that elementary rotation matrices are generated as expected and also work
        as expected.
        """
        rotAngle = 0.0
        rotMatrix = elementary_rotation_matrix("x", rotAngle)
        np.testing.assert_array_almost_equal(np.array([1.0, 0.0, 0.0]), rotMatrix[0, :])
        np.testing.assert_array_almost_equal(np.array([0.0, 1.0, 0.0]), rotMatrix[1, :])
        np.testing.assert_array_almost_equal(np.array([0.0, 0.0, 1.0]), rotMatrix[2, :])
        rotAngle = np.pi / 6.0
        rotMatrix = elementary_rotation_matrix("x", rotAngle)
        np.testing.assert_array_almost_equal(np.array([1.0, 0.0, 0.0]), rotMatrix[0, :])
        np.testing.assert_array_almost_equal(
            np.array([0.0, 0.5 * np.sqrt(3.0), 0.5]), rotMatrix[1, :]
        )
        np.testing.assert_array_almost_equal(
            np.array([0.0, -0.5, 0.5 * np.sqrt(3.0)]), rotMatrix[2, :]
        )
        rotAngles = 2.0 * np.pi * rng.uniform(size=10)
        for rotAngle in rotAngles:
            rotMatrix = elementary_rotation_matrix("x", rotAngle)
            np.testing.assert_array_almost_equal(
                np.array([1.0, 0.0, 0.0]), rotMatrix[0, :]
            )
            np.testing.assert_array_almost_equal(
                np.array([0.0, np.cos(rotAngle), np.sin(rotAngle)]), rotMatrix[1, :]
            )
            np.testing.assert_array_almost_equal(
                np.array([0.0, -np.sin(rotAngle), np.cos(rotAngle)]), rotMatrix[2, :]
            )

        rotAngle = 0.0
        rotMatrix = elementary_rotation_matrix("y", rotAngle)
        np.testing.assert_array_almost_equal(np.array([1.0, 0.0, 0.0]), rotMatrix[0, :])
        np.testing.assert_array_almost_equal(np.array([0.0, 1.0, 0.0]), rotMatrix[1, :])
        np.testing.assert_array_almost_equal(np.array([0.0, 0.0, 1.0]), rotMatrix[2, :])
        rotAngle = np.pi / 6.0
        rotMatrix = elementary_rotation_matrix("y", rotAngle)
        np.testing.assert_array_almost_equal(
            np.array([0.5 * np.sqrt(3.0), 0.0, -0.5]), rotMatrix[0, :]
        )
        np.testing.assert_array_almost_equal(np.array([0.0, 1.0, 0.0]), rotMatrix[1, :])
        np.testing.assert_array_almost_equal(
            np.array([0.5, 0.0, 0.5 * np.sqrt(3.0)]), rotMatrix[2, :]
        )
        rotAngles = 2.0 * np.pi * rng.uniform(size=10)
        for rotAngle in rotAngles:
            rotMatrix = elementary_rotation_matrix("y", rotAngle)
            np.testing.assert_array_almost_equal(
                np.array([np.cos(rotAngle), 0.0, -np.sin(rotAngle)]), rotMatrix[0, :]
            )
            np.testing.assert_array_almost_equal(
                np.array(
                    [
                        0.0,
                        1.0,
                        0.0,
                    ]
                ),
                rotMatrix[1, :],
            )
            np.testing.assert_array_almost_equal(
                np.array([np.sin(rotAngle), 0.0, np.cos(rotAngle)]), rotMatrix[2, :]
            )

        rotAngle = 0.0
        rotMatrix = elementary_rotation_matrix("z", rotAngle)
        np.testing.assert_array_almost_equal(np.array([1.0, 0.0, 0.0]), rotMatrix[0, :])
        np.testing.assert_array_almost_equal(np.array([0.0, 1.0, 0.0]), rotMatrix[1, :])
        np.testing.assert_array_almost_equal(np.array([0.0, 0.0, 1.0]), rotMatrix[2, :])
        rotAngle = np.pi / 6.0
        rotMatrix = elementary_rotation_matrix("z", rotAngle)
        np.testing.assert_array_almost_equal(
            np.array([0.5 * np.sqrt(3.0), 0.5, 0.0]), rotMatrix[0, :]
        )
        np.testing.assert_array_almost_equal(
            np.array([-0.5, 0.5 * np.sqrt(3.0), 0.0]), rotMatrix[1, :]
        )
        np.testing.assert_array_almost_equal(np.array([0.0, 0.0, 1.0]), rotMatrix[2, :])
        rotAngles = 2.0 * np.pi * rng.uniform(size=10)
        for rotAngle in rotAngles:
            rotMatrix = elementary_rotation_matrix("z", rotAngle)
            np.testing.assert_array_almost_equal(
                np.array([np.cos(rotAngle), np.sin(rotAngle), 0.0]), rotMatrix[0, :]
            )
            np.testing.assert_array_almost_equal(
                np.array([-np.sin(rotAngle), np.cos(rotAngle), 0.0]), rotMatrix[1, :]
            )
            np.testing.assert_array_almost_equal(
                np.array(
                    [
                        0.0,
                        0.0,
                        1.0,
                    ]
                ),
                rotMatrix[2, :],
            )

        np.testing.assert_raises(Exception, elementary_rotation_matrix, "blah", 0.0)

    def testNormalTriad(self):
        """
        Check the correct working of the normalTriad() function.
        """
        phi = np.array([0.0, np.pi / 2, np.pi, 3 * np.pi / 2, 0.0, 0.0])
        theta = np.array([0.0, 0.0, 0.0, 0.0, np.pi / 2, -np.pi / 2])
        pExpected = np.array(
            [[0, -1, 0, 1, 0, 0], [1, 0, -1, 0, 1, 1], [0, 0, 0, 0, 0, 0]]
        )
        qExpected = np.array(
            [[0, 0, 0, 0, -1, 1], [0, 0, 0, 0, 0, 0], [1, 1, 1, 1, 0, 0]]
        )
        rExpected = np.array(
            [[1, 0, -1, 0, 0, 0], [0, 1, 0, -1, 0, 0], [0, 0, 0, 0, 1, -1]]
        )
        p, q, r = normal_triad(phi, theta)
        np.testing.assert_array_almost_equal(pExpected, p, decimal=15)
        np.testing.assert_array_almost_equal(qExpected, q, decimal=15)
        np.testing.assert_array_almost_equal(rExpected, r, decimal=15)

        phiRandom = 2.0 * np.pi * rng.uniform(size=10)
        thetaRandom = -np.pi / 2.0 + np.pi * rng.uniform(size=10)
        z = np.array([0, 0, 1])
        for phi in phiRandom:
            for theta in thetaRandom:
                p, q, r = normal_triad(phi, theta)
                rExpected = np.array(
                    [
                        np.cos(phi) * np.cos(theta),
                        np.sin(phi) * np.cos(theta),
                        np.sin(theta),
                    ]
                )
                pExpected = np.cross(z, rExpected)
                pExpected = pExpected / np.sqrt(np.dot(pExpected, pExpected))
                qExpected = np.cross(rExpected, pExpected)
                np.testing.assert_array_almost_equal(pExpected, p, decimal=8)
                np.testing.assert_array_almost_equal(qExpected, q, decimal=8)
                np.testing.assert_array_almost_equal(rExpected, r, decimal=8)

    def testPhaseSpaceToAstrometry(self):
        """
        Verify that the phaseSpaceToAstrometry() function works correctly.
        """
        auKmYearPerSec = 4.74047
        vx = np.array([1, 0, -1, 0, 0, 0])
        vy = np.array([0, 1, 0, -1, 0, 0])
        vz = np.array([0, 0, 0, 0, 1, -1])

        x = np.ones_like(vx)
        y = np.zeros_like(vy)
        z = np.zeros_like(vz)
        phiExpected = np.zeros_like(x)
        thetaExpected = np.zeros_like(x)
        parallaxExpected = 1000 * np.ones_like(x)
        muphistarExpected = (
            parallaxExpected / auKmYearPerSec * np.array([0, 1, 0, -1, 0, 0])
        )
        muthetaExpected = (
            parallaxExpected / auKmYearPerSec * np.array([0, 0, 0, 0, 1, -1])
        )
        vradExpected = np.array([1, 0, -1, 0, 0, 0])
        phi, theta, parallax, muphistar, mutheta, vrad = phase_space_to_astrometry(
            x, y, z, vx, vy, vz
        )
        np.testing.assert_array_almost_equal(phiExpected, phi, decimal=3)
        np.testing.assert_array_almost_equal(thetaExpected, theta, decimal=3)
        np.testing.assert_array_almost_equal(parallaxExpected, parallax, decimal=3)
        np.testing.assert_array_almost_equal(muphistarExpected, muphistar, decimal=3)
        np.testing.assert_array_almost_equal(muthetaExpected, mutheta, decimal=3)
        np.testing.assert_array_almost_equal(vradExpected, vrad, decimal=3)

        x = np.zeros_like(vx)
        y = 2.0 * np.ones_like(vy)
        z = np.zeros_like(vz)
        phiExpected = np.zeros_like(x) + np.pi / 2.0
        thetaExpected = np.zeros_like(x)
        parallaxExpected = 500 * np.ones_like(x)
        muphistarExpected = (
            parallaxExpected / auKmYearPerSec * np.array([-1, 0, 1, 0, 0, 0])
        )
        muthetaExpected = (
            parallaxExpected / auKmYearPerSec * np.array([0, 0, 0, 0, 1, -1])
        )
        vradExpected = np.array([0, 1, 0, -1, 0, 0])
        phi, theta, parallax, muphistar, mutheta, vrad = phase_space_to_astrometry(
            x, y, z, vx, vy, vz
        )
        np.testing.assert_array_almost_equal(phiExpected, phi, decimal=3)
        np.testing.assert_array_almost_equal(thetaExpected, theta, decimal=3)
        np.testing.assert_array_almost_equal(parallaxExpected, parallax, decimal=3)
        np.testing.assert_array_almost_equal(muphistarExpected, muphistar, decimal=3)
        np.testing.assert_array_almost_equal(muthetaExpected, mutheta, decimal=3)
        np.testing.assert_array_almost_equal(vradExpected, vrad, decimal=3)

        x = -0.5 * np.ones_like(vx)
        y = np.zeros_like(vy)
        z = np.zeros_like(vz)
        phiExpected = np.zeros_like(x) + np.pi
        thetaExpected = np.zeros_like(x)
        parallaxExpected = 2000 * np.ones_like(x)
        muphistarExpected = (
            parallaxExpected / auKmYearPerSec * np.array([0, -1, 0, 1, 0, 0])
        )
        muthetaExpected = (
            parallaxExpected / auKmYearPerSec * np.array([0, 0, 0, 0, 1, -1])
        )
        vradExpected = np.array([-1, 0, 1, 0, 0, 0])
        phi, theta, parallax, muphistar, mutheta, vrad = phase_space_to_astrometry(
            x, y, z, vx, vy, vz
        )
        np.testing.assert_array_almost_equal(phiExpected, phi, decimal=3)
        np.testing.assert_array_almost_equal(thetaExpected, theta, decimal=3)
        np.testing.assert_array_almost_equal(parallaxExpected, parallax, decimal=3)
        np.testing.assert_array_almost_equal(muphistarExpected, muphistar, decimal=3)
        np.testing.assert_array_almost_equal(muthetaExpected, mutheta, decimal=3)
        np.testing.assert_array_almost_equal(vradExpected, vrad, decimal=3)

        x = np.zeros_like(vx)
        y = -1 * np.ones_like(vy)
        z = np.zeros_like(vz)
        phiExpected = np.zeros_like(x) + 3 * np.pi / 2.0
        thetaExpected = np.zeros_like(x)
        parallaxExpected = 1000 * np.ones_like(x)
        muphistarExpected = (
            parallaxExpected / auKmYearPerSec * np.array([1, 0, -1, 0, 0, 0])
        )
        muthetaExpected = (
            parallaxExpected / auKmYearPerSec * np.array([0, 0, 0, 0, 1, -1])
        )
        vradExpected = np.array([0, -1, 0, 1, 0, 0])
        phi, theta, parallax, muphistar, mutheta, vrad = phase_space_to_astrometry(
            x, y, z, vx, vy, vz
        )
        np.testing.assert_array_almost_equal(phiExpected, phi, decimal=3)
        np.testing.assert_array_almost_equal(thetaExpected, theta, decimal=3)
        np.testing.assert_array_almost_equal(parallaxExpected, parallax, decimal=3)
        np.testing.assert_array_almost_equal(muphistarExpected, muphistar, decimal=3)
        np.testing.assert_array_almost_equal(muthetaExpected, mutheta, decimal=3)
        np.testing.assert_array_almost_equal(vradExpected, vrad, decimal=3)

        x = np.zeros_like(vx)
        y = np.zeros_like(vy)
        z = np.ones_like(vz)
        phiExpected = np.zeros_like(x)
        thetaExpected = np.zeros_like(x) + np.pi / 2.0
        parallaxExpected = 1000 * np.ones_like(x)
        muphistarExpected = (
            parallaxExpected / auKmYearPerSec * np.array([0, 1, 0, -1, 0, 0])
        )
        muthetaExpected = (
            parallaxExpected / auKmYearPerSec * np.array([-1, 0, 1, 0, 0, 0])
        )
        vradExpected = np.array([0, 0, 0, 0, 1, -1])
        phi, theta, parallax, muphistar, mutheta, vrad = phase_space_to_astrometry(
            x, y, z, vx, vy, vz
        )
        np.testing.assert_array_almost_equal(phiExpected, phi, decimal=3)
        np.testing.assert_array_almost_equal(thetaExpected, theta, decimal=3)
        np.testing.assert_array_almost_equal(parallaxExpected, parallax, decimal=3)
        np.testing.assert_array_almost_equal(muphistarExpected, muphistar, decimal=3)
        np.testing.assert_array_almost_equal(muthetaExpected, mutheta, decimal=3)
        np.testing.assert_array_almost_equal(vradExpected, vrad, decimal=3)

        x = np.zeros_like(vx)
        y = np.zeros_like(vy)
        z = -1 * np.ones_like(vz)
        phiExpected = np.zeros_like(x)
        thetaExpected = np.zeros_like(x) - np.pi / 2.0
        parallaxExpected = 1000 * np.ones_like(x)
        muphistarExpected = (
            parallaxExpected / auKmYearPerSec * np.array([0, 1, 0, -1, 0, 0])
        )
        muthetaExpected = (
            parallaxExpected / auKmYearPerSec * np.array([1, 0, -1, 0, 0, 0])
        )
        vradExpected = np.array([0, 0, 0, 0, -1, 1])
        phi, theta, parallax, muphistar, mutheta, vrad = phase_space_to_astrometry(
            x, y, z, vx, vy, vz
        )
        np.testing.assert_array_almost_equal(phiExpected, phi, decimal=3)
        np.testing.assert_array_almost_equal(thetaExpected, theta, decimal=3)
        np.testing.assert_array_almost_equal(parallaxExpected, parallax, decimal=3)
        np.testing.assert_array_almost_equal(muphistarExpected, muphistar, decimal=3)
        np.testing.assert_array_almost_equal(muthetaExpected, mutheta, decimal=3)
        np.testing.assert_array_almost_equal(vradExpected, vrad, decimal=3)

    def testAstrometryToPhaseSpace(self):
        """
        Verify that the astrometryToPhaseSpace() function works correctly.
        """
        auKmYearPerSec = 4.74047
        muphistar = np.array([1, -1, 0, 0, 0, 0]) / auKmYearPerSec
        mutheta = np.array([0, 0, 1, -1, 0, 0]) / auKmYearPerSec
        vrad = np.array([0, 0, 0, 0, 1, -1])

        parallax = np.zeros_like(vrad) + 1000.0
        phi = np.zeros_like(vrad)
        theta = np.zeros_like(vrad)
        xExpected = np.ones_like(vrad)
        yExpected = np.zeros_like(vrad)
        zExpected = np.zeros_like(vrad)
        vxExpected = np.array([0, 0, 0, 0, 1, -1])
        vyExpected = np.array([1, -1, 0, 0, 0, 0]) / parallax
        vzExpected = np.array([0, 0, 1, -1, 0, 0]) / parallax
        x, y, z, vx, vy, vz = astrometry_to_phase_space(
            phi, theta, parallax, muphistar, mutheta, vrad
        )
        np.testing.assert_array_almost_equal(xExpected, x, decimal=6)
        np.testing.assert_array_almost_equal(yExpected, y, decimal=6)
        np.testing.assert_array_almost_equal(zExpected, z, decimal=6)
        np.testing.assert_array_almost_equal(vxExpected, vx, decimal=6)
        np.testing.assert_array_almost_equal(vyExpected, vy, decimal=6)
        np.testing.assert_array_almost_equal(vzExpected, vz, decimal=6)

        parallax = np.zeros_like(vrad) + 100.0
        phi = np.zeros_like(vrad) + np.pi / 2.0
        theta = np.zeros_like(vrad)
        xExpected = np.zeros_like(vrad)
        yExpected = 10.0 * np.ones_like(vrad)
        zExpected = np.zeros_like(vrad)
        vxExpected = np.array([-1, 1, 0, 0, 0, 0]) / parallax
        vyExpected = np.array([0, 0, 0, 0, 1, -1])
        vzExpected = np.array([0, 0, 1, -1, 0, 0]) / parallax
        x, y, z, vx, vy, vz = astrometry_to_phase_space(
            phi, theta, parallax, muphistar, mutheta, vrad
        )
        np.testing.assert_array_almost_equal(xExpected, x, decimal=6)
        np.testing.assert_array_almost_equal(yExpected, y, decimal=6)
        np.testing.assert_array_almost_equal(zExpected, z, decimal=6)
        np.testing.assert_array_almost_equal(vxExpected, vx, decimal=6)
        np.testing.assert_array_almost_equal(vyExpected, vy, decimal=6)
        np.testing.assert_array_almost_equal(vzExpected, vz, decimal=6)

        parallax = np.zeros_like(vrad) + 20.0
        phi = np.zeros_like(vrad) + np.pi
        theta = np.zeros_like(vrad)
        xExpected = -50.0 * np.ones_like(vrad)
        yExpected = np.zeros_like(vrad)
        zExpected = np.zeros_like(vrad)
        vxExpected = np.array([0, 0, 0, 0, -1, 1])
        vyExpected = np.array([-1, 1, 0, 0, 0, 0]) / parallax
        vzExpected = np.array([0, 0, 1, -1, 0, 0]) / parallax
        x, y, z, vx, vy, vz = astrometry_to_phase_space(
            phi, theta, parallax, muphistar, mutheta, vrad
        )
        np.testing.assert_array_almost_equal(xExpected, x, decimal=6)
        np.testing.assert_array_almost_equal(yExpected, y, decimal=6)
        np.testing.assert_array_almost_equal(zExpected, z, decimal=6)
        np.testing.assert_array_almost_equal(vxExpected, vx, decimal=6)
        np.testing.assert_array_almost_equal(vyExpected, vy, decimal=6)
        np.testing.assert_array_almost_equal(vzExpected, vz, decimal=6)

        parallax = np.zeros_like(vrad) + 20.0
        phi = np.zeros_like(vrad) - np.pi / 2.0
        theta = np.zeros_like(vrad)
        xExpected = np.zeros_like(vrad)
        yExpected = -50.0 * np.ones_like(vrad)
        zExpected = np.zeros_like(vrad)
        vxExpected = np.array([1, -1, 0, 0, 0, 0]) / parallax
        vyExpected = np.array([0, 0, 0, 0, -1, 1])
        vzExpected = np.array([0, 0, 1, -1, 0, 0]) / parallax
        x, y, z, vx, vy, vz = astrometry_to_phase_space(
            phi, theta, parallax, muphistar, mutheta, vrad
        )
        np.testing.assert_array_almost_equal(xExpected, x, decimal=6)
        np.testing.assert_array_almost_equal(yExpected, y, decimal=6)
        np.testing.assert_array_almost_equal(zExpected, z, decimal=6)
        np.testing.assert_array_almost_equal(vxExpected, vx, decimal=6)
        np.testing.assert_array_almost_equal(vyExpected, vy, decimal=6)
        np.testing.assert_array_almost_equal(vzExpected, vz, decimal=6)

        parallax = np.zeros_like(vrad) + 20.0
        phi = np.zeros_like(vrad)
        theta = np.zeros_like(vrad) + np.pi / 2.0
        xExpected = np.zeros_like(vrad)
        yExpected = np.zeros_like(vrad)
        zExpected = 50.0 * np.ones_like(vrad)
        vxExpected = np.array([0, 0, -1, 1, 0, 0]) / parallax
        vyExpected = np.array([1, -1, 0, 0, 0, 0]) / parallax
        vzExpected = np.array([0, 0, 0, 0, 1, -1])
        x, y, z, vx, vy, vz = astrometry_to_phase_space(
            phi, theta, parallax, muphistar, mutheta, vrad
        )
        np.testing.assert_array_almost_equal(xExpected, x, decimal=6)
        np.testing.assert_array_almost_equal(yExpected, y, decimal=6)
        np.testing.assert_array_almost_equal(zExpected, z, decimal=6)
        np.testing.assert_array_almost_equal(vxExpected, vx, decimal=6)
        np.testing.assert_array_almost_equal(vyExpected, vy, decimal=6)
        np.testing.assert_array_almost_equal(vzExpected, vz, decimal=6)

        parallax = np.zeros_like(vrad) + 20.0
        phi = np.zeros_like(vrad)
        theta = np.zeros_like(vrad) - np.pi / 2.0
        xExpected = np.zeros_like(vrad)
        yExpected = np.zeros_like(vrad)
        zExpected = -50.0 * np.ones_like(vrad)
        vxExpected = np.array([0, 0, 1, -1, 0, 0]) / parallax
        vyExpected = np.array([1, -1, 0, 0, 0, 0]) / parallax
        vzExpected = np.array([0, 0, 0, 0, -1, 1])
        x, y, z, vx, vy, vz = astrometry_to_phase_space(
            phi, theta, parallax, muphistar, mutheta, vrad
        )
        np.testing.assert_array_almost_equal(xExpected, x, decimal=6)
        np.testing.assert_array_almost_equal(yExpected, y, decimal=6)
        np.testing.assert_array_almost_equal(zExpected, z, decimal=6)
        np.testing.assert_array_almost_equal(vxExpected, vx, decimal=6)
        np.testing.assert_array_almost_equal(vyExpected, vy, decimal=6)
        np.testing.assert_array_almost_equal(vzExpected, vz, decimal=6)

        parallax = np.zeros_like(vrad) + 20.0
        phi = np.zeros_like(vrad) + np.pi / 2.0
        theta = np.zeros_like(vrad) - np.pi / 2.0
        xExpected = np.zeros_like(vrad)
        yExpected = np.zeros_like(vrad)
        zExpected = -50.0 * np.ones_like(vrad)
        vxExpected = np.array([-1, 1, 0, 0, 0, 0]) / parallax
        vyExpected = np.array([0, 0, 1, -1, 0, 0]) / parallax
        vzExpected = np.array([0, 0, 0, 0, -1, 1])
        x, y, z, vx, vy, vz = astrometry_to_phase_space(
            phi, theta, parallax, muphistar, mutheta, vrad
        )
        np.testing.assert_array_almost_equal(xExpected, x, decimal=6)
        np.testing.assert_array_almost_equal(yExpected, y, decimal=6)
        np.testing.assert_array_almost_equal(zExpected, z, decimal=6)
        np.testing.assert_array_almost_equal(vxExpected, vx, decimal=6)
        np.testing.assert_array_almost_equal(vyExpected, vy, decimal=6)
        np.testing.assert_array_almost_equal(vzExpected, vz, decimal=6)
