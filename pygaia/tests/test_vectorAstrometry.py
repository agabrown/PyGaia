"""
Unit tests for the vectorastrometry module.
"""

from numpy.testing import TestCase, assert_almost_equal, assert_array_almost_equal, assert_raises
from numpy import pi, array, sqrt, arcsin, cos, sin, cross, dot, zeros_like, ones_like
from numpy.random import rand

from pygaia.astrometry.vectorastrometry import cartesian_to_spherical
from pygaia.astrometry.vectorastrometry import spherical_to_cartesian
from pygaia.astrometry.vectorastrometry import elementary_rotation_matrix
from pygaia.astrometry.vectorastrometry import normal_triad
from pygaia.astrometry.vectorastrometry import phase_space_to_astrometry
from pygaia.astrometry.vectorastrometry import astrometry_to_phase_space

class test_vectorAstrometry(TestCase):

  def test_cartesianToSpherical(self):
    """
    Check correct working of the Cartesian to Spherical coordinates transformation.
    """
    r, phi, theta = cartesian_to_spherical(array([1, 0, 0]), array([0, 1, 0]), array([0, 0, 1]))
    assert_array_almost_equal(r,array([1.0,1.0,1.0]))
    assert_array_almost_equal(phi,array([0.0,pi/2.0,0.0]))
    assert_array_almost_equal(theta,array([0.0,0.0,pi/2.0]))

    r, phi, theta = cartesian_to_spherical(1.0, 1.0, 0.0)
    assert_almost_equal(r, sqrt(2.0))
    assert_almost_equal(phi, pi/4.0)
    assert_almost_equal(theta, 0.0)

    r, phi, theta = cartesian_to_spherical(1.0, 1.0, 1.0)
    assert_almost_equal(r, sqrt(3.0))
    assert_almost_equal(phi, pi/4.0)
    assert_almost_equal(theta, arcsin(1.0/sqrt(3.0)))

    r, phi, theta = cartesian_to_spherical(1.0, 1.0, -1.0)
    assert_almost_equal(r, sqrt(3.0))
    assert_almost_equal(phi, pi/4.0)
    assert_almost_equal(theta, -arcsin(1.0/sqrt(3.0)))

    r, phi, theta = cartesian_to_spherical(-1.0, -1.0, -1.0)
    assert_almost_equal(r, sqrt(3.0))
    assert_almost_equal(phi, pi/4.0+pi)
    assert_almost_equal(theta, -arcsin(1.0/sqrt(3.0)))

    assert_raises(Exception, cartesian_to_spherical, 0.0, 0.0, 0.0)
    assert_raises(Exception, cartesian_to_spherical, array([1, 0, 0, 0]), array([0, 1, 0, 0]), array([0, 0, 0, 1]))

  def test_sphericalToCartesian(self):
    """
    Check correct working of the Spherical to Cartesian coordinates transformation.
    """
    x, y, z = spherical_to_cartesian(array([1.0, 2.0, 3.0, 5.0]), array([0.0, pi / 2.0, 0.0, 3.0 * pi / 4.0]),
                                     array([0.0, 0.0, pi/2.0, pi/6.0]))
    assert_array_almost_equal(x, array([1.0,0.0,0.0, -5.0*sqrt(6.0)/4.0]))
    assert_array_almost_equal(y, array([0.0,2.0,0.0, 5.0*sqrt(6.0)/4.0]))
    assert_array_almost_equal(z, array([0.0,0.0,3.0, 2.5]))

    x, y, z = spherical_to_cartesian(3.0, -pi / 6.0, -pi / 4.0)
    assert_almost_equal(x, 3.0*sqrt(6.0)/4.0)
    assert_almost_equal(y, -3.0*sqrt(2.0)/4.0)
    assert_almost_equal(z, -3.0*sqrt(2.0)/2.0)

  def test_elementaryRotationMatrix(self):
    """
    Check that elementary rotation matrices are generated as expected and also work as expected.
    """
    rotAngle=0.0
    rotMatrix=elementary_rotation_matrix("x", rotAngle)
    assert_array_almost_equal(array([1.0,0.0,0.0]), rotMatrix[0,:])
    assert_array_almost_equal(array([0.0,1.0,0.0]), rotMatrix[1,:])
    assert_array_almost_equal(array([0.0,0.0,1.0]), rotMatrix[2,:])
    rotAngle=pi/6.0
    rotMatrix=elementary_rotation_matrix("x", rotAngle)
    assert_array_almost_equal(array([1.0,0.0,0.0]), rotMatrix[0,:])
    assert_array_almost_equal(array([0.0,0.5*sqrt(3.0),0.5]), rotMatrix[1,:])
    assert_array_almost_equal(array([0.0,-0.5,0.5*sqrt(3.0)]), rotMatrix[2,:])
    rotAngles=2.0*pi*rand(10)
    for rotAngle in rotAngles:
      rotMatrix=elementary_rotation_matrix("x", rotAngle)
      assert_array_almost_equal(array([1.0,0.0,0.0]), rotMatrix[0,:])
      assert_array_almost_equal(array([0.0,cos(rotAngle),sin(rotAngle)]), rotMatrix[1,:])
      assert_array_almost_equal(array([0.0,-sin(rotAngle),cos(rotAngle)]), rotMatrix[2,:])

    rotAngle=0.0
    rotMatrix=elementary_rotation_matrix("y", rotAngle)
    assert_array_almost_equal(array([1.0,0.0,0.0]), rotMatrix[0,:])
    assert_array_almost_equal(array([0.0,1.0,0.0]), rotMatrix[1,:])
    assert_array_almost_equal(array([0.0,0.0,1.0]), rotMatrix[2,:])
    rotAngle=pi/6.0
    rotMatrix=elementary_rotation_matrix("y", rotAngle)
    assert_array_almost_equal(array([0.5*sqrt(3.0),0.0,-0.5]), rotMatrix[0,:])
    assert_array_almost_equal(array([0.0,1.0,0.0]), rotMatrix[1,:])
    assert_array_almost_equal(array([0.5,0.0,0.5*sqrt(3.0)]), rotMatrix[2,:])
    rotAngles=2.0*pi*rand(10)
    for rotAngle in rotAngles:
      rotMatrix=elementary_rotation_matrix("y", rotAngle)
      assert_array_almost_equal(array([cos(rotAngle),0.0,-sin(rotAngle)]), rotMatrix[0,:])
      assert_array_almost_equal(array([0.0,1.0,0.0,]),rotMatrix[1,:])
      assert_array_almost_equal(array([sin(rotAngle),0.0,cos(rotAngle)]), rotMatrix[2,:])

    rotAngle=0.0
    rotMatrix=elementary_rotation_matrix("z", rotAngle)
    assert_array_almost_equal(array([1.0,0.0,0.0]), rotMatrix[0,:])
    assert_array_almost_equal(array([0.0,1.0,0.0]), rotMatrix[1,:])
    assert_array_almost_equal(array([0.0,0.0,1.0]), rotMatrix[2,:])
    rotAngle=pi/6.0
    rotMatrix=elementary_rotation_matrix("z", rotAngle)
    assert_array_almost_equal(array([0.5*sqrt(3.0),0.5,0.0]), rotMatrix[0,:])
    assert_array_almost_equal(array([-0.5,0.5*sqrt(3.0),0.0]), rotMatrix[1,:])
    assert_array_almost_equal(array([0.0,0.0,1.0]), rotMatrix[2,:])
    rotAngles=2.0*pi*rand(10)
    for rotAngle in rotAngles:
      rotMatrix=elementary_rotation_matrix("z", rotAngle)
      assert_array_almost_equal(array([cos(rotAngle),sin(rotAngle),0.0]), rotMatrix[0,:])
      assert_array_almost_equal(array([-sin(rotAngle),cos(rotAngle),0.0]), rotMatrix[1,:])
      assert_array_almost_equal(array([0.0,0.0,1.0,]),rotMatrix[2,:])

    assert_raises(Exception, elementary_rotation_matrix, "blah", 0.0)

  def testNormalTriad(self):
    """
    Check the correct working of the normalTriad() function.
    """
    phi=array([0.0,pi/2,pi,3*pi/2,0.0,0.0])
    theta=array([0.0,0.0,0.0,0.0,pi/2,-pi/2])
    pExpected=array([[0,-1,0,1,0,0], [1,0,-1,0,1,1], [0,0,0,0,0,0]])
    qExpected=array([[0,0,0,0,-1,1], [0,0,0,0,0,0], [1,1,1,1,0,0]])
    rExpected=array([[1,0,-1,0,0,0], [0,1,0,-1,0,0], [0,0,0,0,1,-1]])
    p, q, r = normal_triad(phi, theta)
    assert_array_almost_equal(pExpected, p, decimal=15)
    assert_array_almost_equal(qExpected, q, decimal=15)
    assert_array_almost_equal(rExpected, r, decimal=15)

    phiRandom = 2.0*pi*rand(10)
    thetaRandom = -pi/2.0+pi*rand(10)
    z=array([0,0,1])
    for phi in phiRandom:
      for theta in thetaRandom:
        p, q, r = normal_triad(phi, theta)
        rExpected=array([cos(phi)*cos(theta), sin(phi)*cos(theta), sin(theta)])
        pExpected=cross(z,rExpected)
        pExpected=pExpected/sqrt(dot(pExpected,pExpected))
        qExpected=cross(rExpected,pExpected)
        assert_array_almost_equal(pExpected, p, decimal=8)
        assert_array_almost_equal(qExpected, q, decimal=8)
        assert_array_almost_equal(rExpected, r, decimal=8)

  def testPhaseSpaceToAstrometry(self):
    """
    Verify that the phaseSpaceToAstrometry() function works correctly.
    """
    auKmYearPerSec = 4.74047
    vx = array([1,0,-1,0,0,0])
    vy = array([0,1,0,-1,0,0])
    vz = array([0,0,0,0,1,-1])

    x = ones_like(vx)
    y = zeros_like(vy)
    z = zeros_like(vz)
    phiExpected = zeros_like(x)
    thetaExpected = zeros_like(x)
    parallaxExpected = 1000*ones_like(x)
    muphistarExpected = parallaxExpected/auKmYearPerSec*array([0,1,0,-1,0,0])
    muthetaExpected = parallaxExpected/auKmYearPerSec*array([0,0,0,0,1,-1])
    vradExpected = array([1,0,-1,0,0,0])
    phi, theta, parallax, muphistar, mutheta, vrad = phase_space_to_astrometry(x, y, z, vx, vy, vz)
    assert_array_almost_equal(phiExpected, phi, decimal=3)
    assert_array_almost_equal(thetaExpected, theta, decimal=3)
    assert_array_almost_equal(parallaxExpected, parallax, decimal=3)
    assert_array_almost_equal(muphistarExpected, muphistar, decimal=3)
    assert_array_almost_equal(muthetaExpected, mutheta, decimal=3)
    assert_array_almost_equal(vradExpected, vrad, decimal=3)

    x = zeros_like(vx)
    y = 2.0*ones_like(vy)
    z = zeros_like(vz)
    phiExpected = zeros_like(x)+pi/2.0
    thetaExpected = zeros_like(x)
    parallaxExpected = 500*ones_like(x)
    muphistarExpected = parallaxExpected/auKmYearPerSec*array([-1,0,1,0,0,0])
    muthetaExpected = parallaxExpected/auKmYearPerSec*array([0,0,0,0,1,-1])
    vradExpected = array([0,1,0,-1,0,0])
    phi, theta, parallax, muphistar, mutheta, vrad = phase_space_to_astrometry(x, y, z, vx, vy, vz)
    assert_array_almost_equal(phiExpected, phi, decimal=3)
    assert_array_almost_equal(thetaExpected, theta, decimal=3)
    assert_array_almost_equal(parallaxExpected, parallax, decimal=3)
    assert_array_almost_equal(muphistarExpected, muphistar, decimal=3)
    assert_array_almost_equal(muthetaExpected, mutheta, decimal=3)
    assert_array_almost_equal(vradExpected, vrad, decimal=3)

    x = -0.5*ones_like(vx)
    y = zeros_like(vy)
    z = zeros_like(vz)
    phiExpected = zeros_like(x)+pi
    thetaExpected = zeros_like(x)
    parallaxExpected = 2000*ones_like(x)
    muphistarExpected = parallaxExpected/auKmYearPerSec*array([0,-1,0,1,0,0])
    muthetaExpected = parallaxExpected/auKmYearPerSec*array([0,0,0,0,1,-1])
    vradExpected = array([-1,0,1,0,0,0])
    phi, theta, parallax, muphistar, mutheta, vrad = phase_space_to_astrometry(x, y, z, vx, vy, vz)
    assert_array_almost_equal(phiExpected, phi, decimal=3)
    assert_array_almost_equal(thetaExpected, theta, decimal=3)
    assert_array_almost_equal(parallaxExpected, parallax, decimal=3)
    assert_array_almost_equal(muphistarExpected, muphistar, decimal=3)
    assert_array_almost_equal(muthetaExpected, mutheta, decimal=3)
    assert_array_almost_equal(vradExpected, vrad, decimal=3)

    x = zeros_like(vx)
    y = -1*ones_like(vy)
    z = zeros_like(vz)
    phiExpected = zeros_like(x)+3*pi/2.0
    thetaExpected = zeros_like(x)
    parallaxExpected = 1000*ones_like(x)
    muphistarExpected = parallaxExpected/auKmYearPerSec*array([1,0,-1,0,0,0])
    muthetaExpected = parallaxExpected/auKmYearPerSec*array([0,0,0,0,1,-1])
    vradExpected = array([0,-1,0,1,0,0])
    phi, theta, parallax, muphistar, mutheta, vrad = phase_space_to_astrometry(x, y, z, vx, vy, vz)
    assert_array_almost_equal(phiExpected, phi, decimal=3)
    assert_array_almost_equal(thetaExpected, theta, decimal=3)
    assert_array_almost_equal(parallaxExpected, parallax, decimal=3)
    assert_array_almost_equal(muphistarExpected, muphistar, decimal=3)
    assert_array_almost_equal(muthetaExpected, mutheta, decimal=3)
    assert_array_almost_equal(vradExpected, vrad, decimal=3)

    x = zeros_like(vx)
    y = zeros_like(vy)
    z = ones_like(vz)
    phiExpected = zeros_like(x)
    thetaExpected = zeros_like(x)+pi/2.0
    parallaxExpected = 1000*ones_like(x)
    muphistarExpected = parallaxExpected/auKmYearPerSec*array([0,1,0,-1,0,0])
    muthetaExpected = parallaxExpected/auKmYearPerSec*array([-1,0,1,0,0,0])
    vradExpected = array([0,0,0,0,1,-1])
    phi, theta, parallax, muphistar, mutheta, vrad = phase_space_to_astrometry(x, y, z, vx, vy, vz)
    assert_array_almost_equal(phiExpected, phi, decimal=3)
    assert_array_almost_equal(thetaExpected, theta, decimal=3)
    assert_array_almost_equal(parallaxExpected, parallax, decimal=3)
    assert_array_almost_equal(muphistarExpected, muphistar, decimal=3)
    assert_array_almost_equal(muthetaExpected, mutheta, decimal=3)
    assert_array_almost_equal(vradExpected, vrad, decimal=3)

    x = zeros_like(vx)
    y = zeros_like(vy)
    z = -1*ones_like(vz)
    phiExpected = zeros_like(x)
    thetaExpected = zeros_like(x)-pi/2.0
    parallaxExpected = 1000*ones_like(x)
    muphistarExpected = parallaxExpected/auKmYearPerSec*array([0,1,0,-1,0,0])
    muthetaExpected = parallaxExpected/auKmYearPerSec*array([1,0,-1,0,0,0])
    vradExpected = array([0,0,0,0,-1,1])
    phi, theta, parallax, muphistar, mutheta, vrad = phase_space_to_astrometry(x, y, z, vx, vy, vz)
    assert_array_almost_equal(phiExpected, phi, decimal=3)
    assert_array_almost_equal(thetaExpected, theta, decimal=3)
    assert_array_almost_equal(parallaxExpected, parallax, decimal=3)
    assert_array_almost_equal(muphistarExpected, muphistar, decimal=3)
    assert_array_almost_equal(muthetaExpected, mutheta, decimal=3)
    assert_array_almost_equal(vradExpected, vrad, decimal=3)

  def testAstrometryToPhaseSpace(self):
    """
    Verify that the astrometryToPhaseSpace() function works correctly.
    """
    auKmYearPerSec = 4.74047
    muphistar = array([1,-1,0,0,0,0])/auKmYearPerSec
    mutheta = array([0,0,1,-1,0,0])/auKmYearPerSec
    vrad = array([0,0,0,0,1,-1])

    parallax = zeros_like(vrad)+1000.0
    phi = zeros_like(vrad)
    theta = zeros_like(vrad)
    xExpected = ones_like(vrad)
    yExpected = zeros_like(vrad)
    zExpected = zeros_like(vrad)
    vxExpected = array([0,0,0,0,1,-1])
    vyExpected = array([1,-1,0,0,0,0])/parallax
    vzExpected = array([0,0,1,-1,0,0])/parallax
    x, y, z, vx, vy, vz = astrometry_to_phase_space(phi, theta, parallax, muphistar, mutheta, vrad)
    assert_array_almost_equal(xExpected, x, decimal=6)
    assert_array_almost_equal(yExpected, y, decimal=6)
    assert_array_almost_equal(zExpected, z, decimal=6)
    assert_array_almost_equal(vxExpected, vx, decimal=6)
    assert_array_almost_equal(vyExpected, vy, decimal=6)
    assert_array_almost_equal(vzExpected, vz, decimal=6)

    parallax = zeros_like(vrad)+100.0
    phi = zeros_like(vrad)+pi/2.0
    theta = zeros_like(vrad)
    xExpected = zeros_like(vrad)
    yExpected = 10.0*ones_like(vrad)
    zExpected = zeros_like(vrad)
    vxExpected = array([-1,1,0,0,0,0])/parallax
    vyExpected = array([0,0,0,0,1,-1])
    vzExpected = array([0,0,1,-1,0,0])/parallax
    x, y, z, vx, vy, vz = astrometry_to_phase_space(phi, theta, parallax, muphistar, mutheta, vrad)
    assert_array_almost_equal(xExpected, x, decimal=6)
    assert_array_almost_equal(yExpected, y, decimal=6)
    assert_array_almost_equal(zExpected, z, decimal=6)
    assert_array_almost_equal(vxExpected, vx, decimal=6)
    assert_array_almost_equal(vyExpected, vy, decimal=6)
    assert_array_almost_equal(vzExpected, vz, decimal=6)

    parallax = zeros_like(vrad)+20.0
    phi = zeros_like(vrad)+pi
    theta = zeros_like(vrad)
    xExpected = -50.0*ones_like(vrad)
    yExpected = zeros_like(vrad)
    zExpected = zeros_like(vrad)
    vxExpected = array([0,0,0,0,-1,1])
    vyExpected = array([-1,1,0,0,0,0])/parallax
    vzExpected = array([0,0,1,-1,0,0])/parallax
    x, y, z, vx, vy, vz = astrometry_to_phase_space(phi, theta, parallax, muphistar, mutheta, vrad)
    assert_array_almost_equal(xExpected, x, decimal=6)
    assert_array_almost_equal(yExpected, y, decimal=6)
    assert_array_almost_equal(zExpected, z, decimal=6)
    assert_array_almost_equal(vxExpected, vx, decimal=6)
    assert_array_almost_equal(vyExpected, vy, decimal=6)
    assert_array_almost_equal(vzExpected, vz, decimal=6)

    parallax = zeros_like(vrad)+20.0
    phi = zeros_like(vrad)-pi/2.0
    theta = zeros_like(vrad)
    xExpected = zeros_like(vrad)
    yExpected = -50.0*ones_like(vrad)
    zExpected = zeros_like(vrad)
    vxExpected = array([1,-1,0,0,0,0])/parallax
    vyExpected = array([0,0,0,0,-1,1])
    vzExpected = array([0,0,1,-1,0,0])/parallax
    x, y, z, vx, vy, vz = astrometry_to_phase_space(phi, theta, parallax, muphistar, mutheta, vrad)
    assert_array_almost_equal(xExpected, x, decimal=6)
    assert_array_almost_equal(yExpected, y, decimal=6)
    assert_array_almost_equal(zExpected, z, decimal=6)
    assert_array_almost_equal(vxExpected, vx, decimal=6)
    assert_array_almost_equal(vyExpected, vy, decimal=6)
    assert_array_almost_equal(vzExpected, vz, decimal=6)

    parallax = zeros_like(vrad)+20.0
    phi = zeros_like(vrad)
    theta = zeros_like(vrad)+pi/2.0
    xExpected = zeros_like(vrad)
    yExpected = zeros_like(vrad)
    zExpected = 50.0*ones_like(vrad)
    vxExpected = array([0,0,-1,1,0,0])/parallax
    vyExpected = array([1,-1,0,0,0,0])/parallax
    vzExpected = array([0,0,0,0,1,-1])
    x, y, z, vx, vy, vz = astrometry_to_phase_space(phi, theta, parallax, muphistar, mutheta, vrad)
    assert_array_almost_equal(xExpected, x, decimal=6)
    assert_array_almost_equal(yExpected, y, decimal=6)
    assert_array_almost_equal(zExpected, z, decimal=6)
    assert_array_almost_equal(vxExpected, vx, decimal=6)
    assert_array_almost_equal(vyExpected, vy, decimal=6)
    assert_array_almost_equal(vzExpected, vz, decimal=6)

    parallax = zeros_like(vrad)+20.0
    phi = zeros_like(vrad)
    theta = zeros_like(vrad)-pi/2.0
    xExpected = zeros_like(vrad)
    yExpected = zeros_like(vrad)
    zExpected = -50.0*ones_like(vrad)
    vxExpected = array([0,0,1,-1,0,0])/parallax
    vyExpected = array([1,-1,0,0,0,0])/parallax
    vzExpected = array([0,0,0,0,-1,1])
    x, y, z, vx, vy, vz = astrometry_to_phase_space(phi, theta, parallax, muphistar, mutheta, vrad)
    assert_array_almost_equal(xExpected, x, decimal=6)
    assert_array_almost_equal(yExpected, y, decimal=6)
    assert_array_almost_equal(zExpected, z, decimal=6)
    assert_array_almost_equal(vxExpected, vx, decimal=6)
    assert_array_almost_equal(vyExpected, vy, decimal=6)
    assert_array_almost_equal(vzExpected, vz, decimal=6)

    parallax = zeros_like(vrad)+20.0
    phi = zeros_like(vrad)+pi/2.0
    theta = zeros_like(vrad)-pi/2.0
    xExpected = zeros_like(vrad)
    yExpected = zeros_like(vrad)
    zExpected = -50.0*ones_like(vrad)
    vxExpected = array([-1,1,0,0,0,0])/parallax
    vyExpected = array([0,0,1,-1,0,0])/parallax
    vzExpected = array([0,0,0,0,-1,1])
    x, y, z, vx, vy, vz = astrometry_to_phase_space(phi, theta, parallax, muphistar, mutheta, vrad)
    assert_array_almost_equal(xExpected, x, decimal=6)
    assert_array_almost_equal(yExpected, y, decimal=6)
    assert_array_almost_equal(zExpected, z, decimal=6)
    assert_array_almost_equal(vxExpected, vx, decimal=6)
    assert_array_almost_equal(vyExpected, vy, decimal=6)
    assert_array_almost_equal(vzExpected, vz, decimal=6)
