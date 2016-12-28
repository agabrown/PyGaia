"""
Unit tests for the errors.astrometric module.
"""

from numpy.testing import TestCase, assert_almost_equal, assert_array_almost_equal, assert_raises
from numpy import pi, linspace, power

from pygaia.errors import astrometric as astrom

class test_errorsAstrometric(TestCase):

  def test_errorScalingFactor(self):
    """
    Check that this function works for arrays and scalars.
    """
    beta = linspace(0,pi/2.0,30)
    factors = astrom.errorScalingFactor("parallax", beta)
    for factorValue in factors:
      self.assertTrue(factorValue>0.0)
    self.assertTrue(astrom.errorScalingFactor("parallax", pi/3.0)>0.0)

  def test_parallaxErrorSkyAvg(self):
    """
    Verify that the function works for reasonable combinations of array and scalar input parameters.
    """
    self.assertTrue(astrom.parallaxErrorSkyAvg(15.0,3.0)>0.0)
    gmags = linspace(6,20,100)
    vmini = linspace(-1,4,100)
    errors = astrom.parallaxErrorSkyAvg(gmags, vmini)
    for error in errors:
      self.assertTrue(error>0.0)
    errors = astrom.parallaxErrorSkyAvg(gmags, 1.0)
    for error in errors:
      self.assertTrue(error>0.0)
    errors = astrom.parallaxErrorSkyAvg(20.0, vmini)
    for error in errors:
      self.assertTrue(error>0.0)

  def test_parallaxErrorSkyAveExtended(self):
      """
      Check that errors are correctly scaled for a mission extension.
      """
      self.assertTrue(astrom.parallaxErrorSkyAvg(15.0,3.0)>astrom.parallaxErrorSkyAvg(15.0,3.0,extension=1.0))
      error = astrom.parallaxErrorSkyAvg(15.0,3.0)
      for e in range(6):
          exterror = astrom.parallaxErrorSkyAvg(15.0,3.0,extension=e)
          assert_almost_equal(exterror/error, power((5.0+e)/5.0,-0.5), decimal=8)

  def test_parallaxError(self):
    """
    Verify that the function works for reasonable combinations of array and scalar input parameters.
    """
    self.assertTrue(astrom.parallaxError(15.0,3.0,pi/2.0)>0.0)
    gmags = linspace(6,20,100)
    vmini = linspace(-1,4,100)
    beta = linspace(0,pi/2.0,100)
    errors = astrom.parallaxError(gmags, vmini, beta)
    for error in errors:
      self.assertTrue(error>0.0)
    errors = astrom.parallaxError(gmags, 1.0, beta)
    for error in errors:
      self.assertTrue(error>0.0)
    errors = astrom.parallaxError(20.0, vmini, beta)
    for error in errors:
      self.assertTrue(error>0.0)
    errors = astrom.parallaxError(gmags, vmini, 0.0)
    for error in errors:
      self.assertTrue(error>0.0)

  def test_positionErrorSkyAvg(self):
    """
    Verify that the function works for reasonable combinations of array and scalar input parameters.
    """
    sigalphastar, sigdelta = astrom.positionErrorSkyAvg(15.0,3.0)
    self.assertTrue(sigalphastar>0.0)
    self.assertTrue(sigdelta>0.0)

    gmags = linspace(6,20,100)
    vmini = linspace(-1,4,100)
    sigalphastar, sigdelta = astrom.positionErrorSkyAvg(gmags, vmini)
    for error in sigalphastar:
      self.assertTrue(error>0.0)
    for error in sigdelta:
      self.assertTrue(error>0.0)
    sigalphastar, sigdelta = astrom.positionErrorSkyAvg(gmags, -1.0)
    for error in sigalphastar:
      self.assertTrue(error>0.0)
    for error in sigdelta:
      self.assertTrue(error>0.0)
    sigalphastar, sigdelta = astrom.positionErrorSkyAvg(15.0, vmini)
    for error in sigalphastar:
      self.assertTrue(error>0.0)
    for error in sigdelta:
      self.assertTrue(error>0.0)

  def test_positionErrorSkyAveExtended(self):
      """
      Check that errors are correctly scaled for a mission extension.
      """
      sigalphastar, sigdelta = astrom.positionErrorSkyAvg(15.0,3.0)
      for e in range(6):
          sigalphastarext, sigdeltaext = astrom.positionErrorSkyAvg(15.0,3.0,extension=e)
          assert_almost_equal(sigalphastarext/sigalphastar, power((5.0+e)/5.0,-0.5), decimal=8)
          assert_almost_equal(sigdeltaext/sigdelta, power((5.0+e)/5.0,-0.5), decimal=8)

  def test_positionError(self):
    """
    Verify that the function works for reasonable combinations of array and scalar input parameters.
    """
    sigalphastar, sigdelta = astrom.positionError(15.0,3.0,pi/4.0)
    self.assertTrue(sigalphastar>0.0)
    self.assertTrue(sigdelta>0.0)

    gmags = linspace(6,20,100)
    vmini = linspace(-1,4,100)
    beta = linspace(0,pi/2.0,100)
    sigalphastar, sigdelta = astrom.positionError(gmags, vmini, beta)
    for error in sigalphastar:
      self.assertTrue(error>0.0)
    for error in sigdelta:
      self.assertTrue(error>0.0)
    sigalphastar, sigdelta = astrom.positionError(gmags, -1.0, beta)
    for error in sigalphastar:
      self.assertTrue(error>0.0)
    for error in sigdelta:
      self.assertTrue(error>0.0)
    sigalphastar, sigdelta = astrom.positionError(15.0, vmini, beta)
    for error in sigalphastar:
      self.assertTrue(error>0.0)
    for error in sigdelta:
      self.assertTrue(error>0.0)
    sigalphastar, sigdelta = astrom.positionError(gmags, vmini, 0.0)
    for error in sigalphastar:
      self.assertTrue(error>0.0)
    for error in sigdelta:
      self.assertTrue(error>0.0)

  def test_properMotionErrorSkyAvg(self):
    """
    Verify that the function works for reasonable combinations of array and scalar input parameters.
    """
    sigalphastar, sigdelta = astrom.properMotionErrorSkyAvg(15.0,3.0)
    self.assertTrue(sigalphastar>0.0)
    self.assertTrue(sigdelta>0.0)

    gmags = linspace(6,20,100)
    vmini = linspace(-1,4,100)
    sigalphastar, sigdelta = astrom.properMotionErrorSkyAvg(gmags, vmini)
    for error in sigalphastar:
      self.assertTrue(error>0.0)
    for error in sigdelta:
      self.assertTrue(error>0.0)
    sigalphastar, sigdelta = astrom.properMotionErrorSkyAvg(gmags, -1.0)
    for error in sigalphastar:
      self.assertTrue(error>0.0)
    for error in sigdelta:
      self.assertTrue(error>0.0)
    sigalphastar, sigdelta = astrom.properMotionErrorSkyAvg(15.0, vmini)
    for error in sigalphastar:
      self.assertTrue(error>0.0)
    for error in sigdelta:
      self.assertTrue(error>0.0)

  def test_properMotionErrorSkyAveExtended(self):
      """
      Check that errors are correctly scaled for a mission extension.
      """
      sigmualphastar, sigmudelta = astrom.properMotionErrorSkyAvg(15.0,3.0)
      for e in range(6):
          sigmualphastarext, sigmudeltaext = astrom.properMotionErrorSkyAvg(15.0,3.0,extension=e)
          assert_almost_equal(sigmualphastarext/sigmualphastar, power((5.0+e)/5.0,-1.5), decimal=8)
          assert_almost_equal(sigmudeltaext/sigmudelta, power((5.0+e)/5.0,-1.5), decimal=8)

  def test_properMotionError(self):
    """
    Verify that the function works for reasonable combinations of array and scalar input parameters.
    """
    sigalphastar, sigdelta = astrom.properMotionError(15.0,3.0,pi/4.0)
    self.assertTrue(sigalphastar>0.0)
    self.assertTrue(sigdelta>0.0)

    gmags = linspace(6,20,100)
    vmini = linspace(-1,4,100)
    beta = linspace(0,pi/2.0,100)
    sigalphastar, sigdelta = astrom.properMotionError(gmags, vmini, beta)
    for error in sigalphastar:
      self.assertTrue(error>0.0)
    for error in sigdelta:
      self.assertTrue(error>0.0)
    sigalphastar, sigdelta = astrom.properMotionError(gmags, -1.0, beta)
    for error in sigalphastar:
      self.assertTrue(error>0.0)
    for error in sigdelta:
      self.assertTrue(error>0.0)
    sigalphastar, sigdelta = astrom.properMotionError(15.0, vmini, beta)
    for error in sigalphastar:
      self.assertTrue(error>0.0)
    for error in sigdelta:
      self.assertTrue(error>0.0)
    sigalphastar, sigdelta = astrom.properMotionError(gmags, vmini, 0.0)
    for error in sigalphastar:
      self.assertTrue(error>0.0)
    for error in sigdelta:
      self.assertTrue(error>0.0)
