"""
Unit tests for the utils module.
"""

import numpy as np

from pygaia.utils import construct_covariance_matrix
from pygaia.astrometry.constants import au_km_year_per_sec


class test_utils(np.testing.TestCase):
    def test_construct_covariance_matrix_onesource(self):
        """
        Check that the covariance matrix is constructed correctly from the inputs
        provided.
        """

        cvec = np.arange(15) + 1
        parallax = 16
        vrad = 17 * au_km_year_per_sec
        vrad_error = 18 * au_km_year_per_sec

        cmat = construct_covariance_matrix(cvec, parallax, vrad, vrad_error)

        self.assertEquals(6, cmat.shape[0])
        self.assertEquals(6, cmat.shape[1])

        iu = np.triu_indices(6, k=1)

        np.testing.assert_equal(cmat[iu], cmat[iu[1], iu[0]])
        np.testing.assert_equal((np.arange(5) + 1) ** 2, np.diag(cmat)[0:-1])

        expected = np.array([12, 21, 32, 45, 60, 88, 120, 156, 210, 300])
        iuu = np.triu_indices(5, k=1)
        np.testing.assert_equal(expected, cmat[:-1, :-1][iuu])

        for i in range(5):
            np.testing.assert_almost_equal(cmat[i, 5], cmat[i, 2] * 17, 1)

        np.testing.assert_almost_equal(
            cmat[5, 5], 9 * (17**2 + 18**2) + (16 * 18) ** 2, 1
        )

    def test_construct_covariance_matrix_nsources(self):
        """
        Check that the covariance matrix is constructed correctly from the inputs
        provided.
        """

        cvec = np.zeros((2, 15))
        cvec[0] = np.arange(15) + 1
        cvec[1] = -cvec[0]
        parallax = np.array([16, -16])
        vrad = np.array([17, -17]) * au_km_year_per_sec
        vrad_error = np.array([18, -18]) * au_km_year_per_sec

        cmat = construct_covariance_matrix(cvec, parallax, vrad, vrad_error)

        self.assertEquals(2, cmat.shape[0])
        self.assertEquals(6, cmat.shape[1])
        self.assertEquals(6, cmat.shape[2])

        iu = np.triu_indices(6, k=1)

        np.testing.assert_equal(cmat[0, iu[0], iu[1]], cmat[0, iu[1], iu[0]])
        np.testing.assert_equal((np.arange(5) + 1) ** 2, np.diag(cmat[0])[0:-1])

        expected = np.array([12, 21, 32, 45, 60, 88, 120, 156, 210, 300])
        iuu = np.triu_indices(5, k=1)
        np.testing.assert_equal(expected, cmat[0, :-1, :-1][iuu])

        for i in range(5):
            np.testing.assert_almost_equal(cmat[0, i, 5], cmat[0, i, 2] * 17, 1)

        np.testing.assert_almost_equal(
            cmat[0, 5, 5], 9 * (17**2 + 18**2) + (16 * 18) ** 2, 1
        )

        np.testing.assert_equal(cmat[1, iu[0], iu[1]], cmat[1, iu[1], iu[0]])
        np.testing.assert_equal((np.arange(5) + 1) ** 2, np.diag(cmat[1])[0:-1])

        expected = -np.array([12, 21, 32, 45, 60, 88, 120, 156, 210, 300])
        iuu = np.triu_indices(5, k=1)
        np.testing.assert_equal(expected, cmat[1, :-1, :-1][iuu])

        for i in range(5):
            np.testing.assert_almost_equal(cmat[1, i, 5], cmat[1, i, 2] * (-17), 1)

        np.testing.assert_almost_equal(
            cmat[1, 5, 5], 9 * (17**2 + 18**2) + (16 * 18) ** 2, 1
        )
