"""
Unit tests for the coordinates module.
"""

from numpy.testing import TestCase, assert_allclose, assert_array_almost_equal
from numpy import pi, array, transpose, cos, sin
from numpy.random import rand

from pygaia.astrometry.coordinates import CoordinateTransformation
from pygaia.astrometry.coordinates import Transformations
from pygaia.astrometry.vectorastrometry import cartesianToSpherical

class test_coordinates(TestCase):

  def setUp(self):
    self.basisVectorXValues=array([1,0,0])
    self.basisVectorYValues=array([0,1,0])
    self.basisVectorZValues=array([0,0,1])

    self.expectedMatIcrsToGal = array([[-0.05, -0.87, -0.48], [+0.49, -0.44, +0.75], [-0.87, -0.20,
      +0.46]])
    self.expectedMatEclipticToIcrs = array([[1.0, 0.0, 0.0], [0.0, 0.92, -0.40], [0.0, 0.40, 0.92]])
    self.expectedMatGalacticToEcliptic = array([[-0.05, 0.49, -0.87], [-0.99, -0.11, 0.00], [-0.1, 0.86,
      0.50]])

  def test_galacticToIcrs(self):
    """
    Verify correctness of transformations from the Galactic to ICRS coordinate systems.
    """
    ct = CoordinateTransformation(Transformations.GAL2ICRS)
    x, y, z = ct.transformCartesianCoordinates(self.basisVectorXValues, self.basisVectorYValues,
        self.basisVectorZValues)
    assert_array_almost_equal(x, self.expectedMatIcrsToGal[:,0], decimal=2)
    assert_array_almost_equal(y, self.expectedMatIcrsToGal[:,1], decimal=2)
    assert_array_almost_equal(z, self.expectedMatIcrsToGal[:,2], decimal=2)

    r, expectedAlpha, expectedDelta = cartesianToSpherical(self.expectedMatIcrsToGal[:,0],
        self.expectedMatIcrsToGal[:,1], self.expectedMatIcrsToGal[:,2])
    alpha, delta = ct.transformSkyCoordinates(array([0.0,pi/2.0,0.0]), array([0.0, 0.0, pi/2.0]))
    assert_array_almost_equal(expectedAlpha, alpha, decimal=1)
    assert_array_almost_equal(expectedDelta, delta, decimal=1)

    alpha, delta = ct.transformSkyCoordinates(0.0,pi/2.0)
    assert_allclose(alpha/pi*180, (192.9-360.0), atol=1.0)
    assert_allclose(delta/pi*180, 27.1, atol=1.0)

  def test_icrsToGalactic(self):
    """
    Verify correctness of transformations from the ICRS to Galactic coordinate systems.
    """
    ct=CoordinateTransformation(Transformations.ICRS2GAL)
    x, y, z = ct.transformCartesianCoordinates(self.basisVectorXValues, self.basisVectorYValues,
        self.basisVectorZValues)
    assert_array_almost_equal(x, self.expectedMatIcrsToGal[0,:], decimal=2)
    assert_array_almost_equal(y, self.expectedMatIcrsToGal[1,:], decimal=2)
    assert_array_almost_equal(z, self.expectedMatIcrsToGal[2,:], decimal=2)

    r, expectedGalon, expectedGalat = cartesianToSpherical(self.expectedMatIcrsToGal[0,:],
        self.expectedMatIcrsToGal[1,:], self.expectedMatIcrsToGal[2,:])
    galon, galat = ct.transformSkyCoordinates(array([0.0,pi/2.0,0.0]), array([0.0, 0.0, pi/2.0]))
    assert_array_almost_equal(expectedGalon, galon, decimal=1)
    assert_array_almost_equal(expectedGalat, galat, decimal=1)

  def test_icrsToEcliptic(self):
    """
    Verify correctness of transformations from the ICRS to Ecliptic coordinate systems.
    """
    ct=CoordinateTransformation(Transformations.ICRS2ECL)
    x, y, z = ct.transformCartesianCoordinates(self.basisVectorXValues, self.basisVectorYValues,
        self.basisVectorZValues)
    assert_array_almost_equal(x, self.expectedMatEclipticToIcrs[:,0], decimal=2)
    assert_array_almost_equal(y, self.expectedMatEclipticToIcrs[:,1], decimal=2)
    assert_array_almost_equal(z, self.expectedMatEclipticToIcrs[:,2], decimal=2)

    r, expectedLambda, expectedBeta = cartesianToSpherical(self.expectedMatEclipticToIcrs[:,0],
        self.expectedMatEclipticToIcrs[:,1], self.expectedMatEclipticToIcrs[:,2])
    lambdaEcl, betaEcl = ct.transformSkyCoordinates(array([0.0,pi/2.0,0.0]), array([0.0, 0.0, pi/2.0]))
    assert_array_almost_equal(expectedLambda, lambdaEcl, decimal=1)
    assert_array_almost_equal(expectedBeta, betaEcl, decimal=1)

    alpha = 2.0*pi*rand(100)
    delta = -pi/2.0+pi*rand(100)
    lambdaEcl, betaEcl = ct.transformSkyCoordinates(alpha, delta)
    result = 0.9175*sin(delta)-0.3978*sin(alpha)*cos(delta)
    assert_array_almost_equal(result, sin(betaEcl), decimal=2)

  def test_eclipticToIcrs(self):
    """
    Verify correctness of transformations from the Ecliptic to ICRS coordinate systems.
    """
    ct=CoordinateTransformation(Transformations.ECL2ICRS)
    x, y, z = ct.transformCartesianCoordinates(self.basisVectorXValues, self.basisVectorYValues,
        self.basisVectorZValues)
    assert_array_almost_equal(x, self.expectedMatEclipticToIcrs[0,:], decimal=2)
    assert_array_almost_equal(y, self.expectedMatEclipticToIcrs[1,:], decimal=2)
    assert_array_almost_equal(z, self.expectedMatEclipticToIcrs[2,:], decimal=2)

    r, expectedAlpha, expectedDelta = cartesianToSpherical(self.expectedMatEclipticToIcrs[0,:],
        self.expectedMatEclipticToIcrs[1,:], self.expectedMatEclipticToIcrs[2,:])
    alpha, delta = ct.transformSkyCoordinates(array([0.0,pi/2.0,0.0]), array([0.0, 0.0, pi/2.0]))
    assert_array_almost_equal(expectedAlpha, alpha, decimal=1)
    assert_array_almost_equal(expectedDelta, delta, decimal=1)

    lambdaEcl = 2.0*pi*rand(100)
    betaEcl = -pi/2.0+pi*rand(100)
    alpha, delta = ct.transformSkyCoordinates(lambdaEcl, betaEcl)
    result = 0.9175*sin(betaEcl)+0.3978*sin(lambdaEcl)*cos(betaEcl)
    assert_array_almost_equal(result, sin(delta), decimal=2)

  def test_eclipticToGalactic(self):
    """
    Verify correctness of transformations from the Ecliptic to Galactic coordinate systems.
    """
    ct=CoordinateTransformation(Transformations.ECL2GAL)
    x, y, z = ct.transformCartesianCoordinates(self.basisVectorXValues, self.basisVectorYValues,
        self.basisVectorZValues)
    assert_array_almost_equal(x, self.expectedMatGalacticToEcliptic[:,0], decimal=2)
    assert_array_almost_equal(y, self.expectedMatGalacticToEcliptic[:,1], decimal=2)
    assert_array_almost_equal(z, self.expectedMatGalacticToEcliptic[:,2], decimal=2)

    r, expectedGalon, expectedGalat = cartesianToSpherical(self.expectedMatGalacticToEcliptic[:,0],
        self.expectedMatGalacticToEcliptic[:,1], self.expectedMatGalacticToEcliptic[:,2])
    galon, galat = ct.transformSkyCoordinates(array([0.0,pi/2.0,0.0]), array([0.0, 0.0, pi/2.0]))
    assert_array_almost_equal(expectedGalon, galon, decimal=1)
    assert_array_almost_equal(expectedGalat, galat, decimal=1)

  def test_galacticToEcliptic(self):
    """
    Verify correctness of transformations from the Galactic to Ecliptic coordinate systems.
    """
    ct=CoordinateTransformation(Transformations.GAL2ECL)
    x, y, z = ct.transformCartesianCoordinates(self.basisVectorXValues, self.basisVectorYValues,
        self.basisVectorZValues)
    assert_array_almost_equal(x, self.expectedMatGalacticToEcliptic[0,:], decimal=2)
    assert_array_almost_equal(y, self.expectedMatGalacticToEcliptic[1,:], decimal=2)
    assert_array_almost_equal(z, self.expectedMatGalacticToEcliptic[2,:], decimal=2)

    r, expectedLambda, expectedBeta = cartesianToSpherical(self.expectedMatGalacticToEcliptic[0,:],
        self.expectedMatGalacticToEcliptic[1,:], self.expectedMatGalacticToEcliptic[2,:])
    lambdaEcl, betaEcl = ct.transformSkyCoordinates(array([0.0,pi/2.0,0.0]), array([0.0, 0.0, pi/2.0]))
    #
    # Note the abs() in the line below is required to avoid an error when -pi and pi are compared.
    # The better solution of course is to write an assert function that can handle modulo 2*pi cases.
    #
    assert_array_almost_equal(abs(expectedLambda), abs(lambdaEcl), decimal=1)
    assert_array_almost_equal(expectedBeta, betaEcl, decimal=1)

    galon = 2.0*pi*rand(100)
    galat = -pi/2.0+pi*rand(100)
    lambdaEcl, betaEcl = ct.transformSkyCoordinates(galon, galat)
    phi=6.38/180.0*pi
    result = 0.4971*sin(galat)+0.8677*sin(galon-phi)*cos(galat)
    assert_array_almost_equal(result, sin(betaEcl), decimal=2)
