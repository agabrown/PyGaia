"""
Unit tests for the errors.astrometric module.
"""

import numpy as np
from numpy.testing import TestCase

from pygaia.errors import astrometric as astrom


class TestErrorsAstrometric(TestCase):

    def test_parallax_uncertainty(self):
        """
        Verify that the function works for reasonable combinations of array and scalar input parameters.
        """
        self.assertTrue(astrom.parallax_uncertainty(15.0) > 0.0)
        self.assertTrue(astrom.parallax_uncertainty(15.0, release='dr4') > 0.0)
        self.assertTrue(astrom.parallax_uncertainty(15.0, release='dr3') > 0.0)
        self.assertTrue(astrom.parallax_uncertainty(15.0, release='dr5') > 0.0)
        gmags = np.linspace(6, 20, 100)
        errorsdr4 = astrom.parallax_uncertainty(gmags)
        errorsdr3 = astrom.parallax_uncertainty(gmags, release='dr3')
        errorsdr5 = astrom.parallax_uncertainty(gmags, release='dr5')
        for errordr4, errordr3, errordr5 in zip(errorsdr4, errorsdr3, errorsdr5):
            self.assertTrue(errordr4 > 0.0)
            self.assertTrue(errordr3 > errordr4)
            self.assertTrue(errordr5 < errordr4)

    def test_position_uncertainty(self):
        """
        Verify that the function works for reasonable combinations of array and scalar input parameters.
        """
        sigalphastar, sigdelta = astrom.position_uncertainty(15.0)
        self.assertTrue(sigalphastar > 0.0)
        self.assertTrue(sigdelta > 0.0)
        sigalphastar, sigdelta = astrom.position_uncertainty(15.0, release='dr3')
        self.assertTrue(sigalphastar > 0.0)
        self.assertTrue(sigdelta > 0.0)
        sigalphastar, sigdelta = astrom.position_uncertainty(15.0, release='dr5')
        self.assertTrue(sigalphastar > 0.0)
        self.assertTrue(sigdelta > 0.0)

        gmags = np.linspace(6, 20, 100)
        sigalphastardr4, sigdeltadr4 = astrom.position_uncertainty(gmags)
        sigalphastardr3, sigdeltadr3 = astrom.position_uncertainty(gmags, release='dr3')
        sigalphastardr5, sigdeltadr5 = astrom.position_uncertainty(gmags, release='dr5')
        for errordr4, errordr3, errordr5 in zip(sigalphastardr4, sigalphastardr3, sigalphastardr5):
            self.assertTrue(errordr4 > 0.0)
            self.assertTrue(errordr3 > errordr4)
            self.assertTrue(errordr5 < errordr4)
        for errordr4, errordr3, errordr5 in zip(sigdeltadr4, sigdeltadr3, sigdeltadr5):
            self.assertTrue(errordr4 > 0.0)
            self.assertTrue(errordr3 > errordr4)
            self.assertTrue(errordr5 < errordr4)

    def test_proper_motion_uncertainty(self):
        """
        Verify that the function works for reasonable combinations of array and scalar input parameters.
        """
        sigpmra, sigpmdec = astrom.proper_motion_uncertainty(15.0)
        self.assertTrue(sigpmra > 0.0)
        self.assertTrue(sigpmdec > 0.0)
        sigpmra, sigpmdec = astrom.proper_motion_uncertainty(15.0, release='dr3')
        self.assertTrue(sigpmra > 0.0)
        self.assertTrue(sigpmdec > 0.0)
        sigpmra, sigpmdec = astrom.proper_motion_uncertainty(15.0, release='dr5')
        self.assertTrue(sigpmra > 0.0)
        self.assertTrue(sigpmdec > 0.0)

        gmags = np.linspace(6, 20, 100)
        sigpmradr4, sigpmdecdr4 = astrom.proper_motion_uncertainty(gmags)
        sigpmradr3, sigpmdecdr3 = astrom.proper_motion_uncertainty(gmags, release='dr3')
        sigpmradr5, sigpmdecdr5 = astrom.proper_motion_uncertainty(gmags, release='dr5')
        for errordr4, errordr3, errordr5 in zip(sigpmradr4, sigpmradr3, sigpmradr5):
            self.assertTrue(errordr4 > 0.0)
            self.assertTrue(errordr3 > errordr4)
            self.assertTrue(errordr5 < errordr4)
        for errordr4, errordr3, errordr5 in zip(sigpmdecdr4, sigpmdecdr3, sigpmdecdr5):
            self.assertTrue(errordr4 > 0.0)
            self.assertTrue(errordr3 > errordr4)
            self.assertTrue(errordr5 < errordr4)
