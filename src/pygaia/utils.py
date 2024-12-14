"""
Utility functions for PyGaia.
"""

import numpy as np

from pygaia.astrometry.constants import au_km_year_per_sec
from astropy.time import Time

__all__ = ["construct_covariance_matrix", "gaiadr_timespan"]

_supported_releases = ["dr1", "dr2", "dr3", "dr4", "dr5"]
_default_release = "dr4"
_jdref = Time("2010-01-01T00:00:00", scale="tcb").jd
_dpac_times = {
    "start": {
        "obmt": 23292998211919880,
        "tcb": Time("2014-07-25T10:31:25.554960001", scale="tcb"),
    },
    "end_epsl": {
        "obmt": 25749998222167256,
        "tcb": Time("2014-08-22T21:01:25.599970336", scale="tcb"),
    },
    "dr1": {
        "obmt": 59429199262402272,
        "tcb": Time("2015-09-16T16:21:27.121893186", scale="tcb"),
    },
    "dr2": {
        "obmt": 81012099291189247,
        "tcb": Time("2016-05-23T11:36:27.459006034", scale="tcb"),
    },
    "dr3": {
        "obmt": 112969900338538592,
        "tcb": Time("2017-05-28T08:46:28.954612431", scale="tcb"),
    },
    "dr4": {
        "obmt": 196566400491690272,
        "tcb": Time("2020-01-20T22:01:30.250520158", scale="tcb"),
    },
    "dr5": {
        "obmt": 353907908175187136,
        "tcb": Time("2025-01-15T06:30:00", scale="tcb"),
    },
}


def construct_covariance_matrix(cvec, parallax, radial_velocity, radial_velocity_error):
    """
    Take the astrometric parameter standard uncertainties and the uncertainty
    correlations as quoted in the Gaia catalogue and construct the covariance matrix.

    Parameters
    ----------
    cvec : array_like
        Array of shape (15,) (1 source) or (N,15) (N sources) for the astrometric
        parameter standard uncertainties and their correlations, as listed in the Gaia
        catalogue [ra_error, dec_error, parallax_error, pmra_error, pmdec_error,
        ra_dec_corr, ra_parallax_corr, ra_pmra_corr, ra_pmdec_corr, dec_parallax_corr,
        dec_pmra_corr, dec_pmdec_corr, parallax_pmra_corr, parallax_pmdec_corr,
        pmra_pmdec_corr]. Units are (mas^2, mas^2/yr, mas^2/yr^2).
    parallax : array_like
        Source parallax(es) array of shape (N,) (mas).
    radial_velocity : array_like
        Source radial velocity (km/s, does not have to be from Gaia RVS!) as array of
        shape (N,). If the radial velocity is not known it can be set to zero.
    radial_velocity_error : array_like
        Source radial velocity  uncertainty (km/s), array of shape (N,). If the radial
        velocity is not know this can be set to the radial velocity dispersion for the
        population the source was drawn from.

    Returns
    -------
    cmat : float array
        Covariance matrix as a Nx6x6 array.
    """

    if np.ndim(cvec) == 1:
        cmat = np.zeros((1, 6, 6))
        nsources = 1
        cv = np.atleast_2d(cvec)
    else:
        nsources = cvec.shape[0]
        cmat = np.zeros((nsources, 6, 6))
        cv = cvec
    for k in range(nsources):
        np.fill_diagonal(cmat[k], cv[k, 0:5] ** 2)

    iu = np.triu_indices(5, k=1)
    for k in range(10):
        i = iu[0][k]
        j = iu[1][k]
        cmat[:, i, j] = cv[:, i] * cv[:, j] * cv[:, k + 5]
        cmat[:, j, i] = cmat[:, i, j]

    for k in range(nsources):
        cmat[k, 0:5, 5] = (
            cmat[k, 0:5, 2] * np.atleast_1d(radial_velocity)[k] / au_km_year_per_sec
        )
    cmat[:, 5, 0:5] = cmat[:, 0:5, 5]
    cmat[:, 5, 5] = (
        cmat[:, 2, 2]
        * (radial_velocity**2 + radial_velocity_error**2)
        / au_km_year_per_sec**2
        + (parallax * radial_velocity_error / au_km_year_per_sec) ** 2
    )

    return np.squeeze(cmat)


def _check_release(release):
    """
    Check if the release requested by the user is supported.

    Parameters
    ----------
    release : str
        Release for which performance predictions are requested.

    Raises
    ------
    ValueError
        When an invalid string is specified for the release parameter.
    """
    if not (release in _supported_releases):
        raise ValueError("Release must be one of dr1, dr2, dr3, dr4, dr5")


def gaiadr_timespan(release=_default_release, epsl=True, absolute=True):
    """
    Provides the observation time spans for Gaia data releases. These are based on the exact time boundaries in Gaia's on board mission timeline (OBMT), see `Gaia Collaboration (2016) <https://ui.adsabs.harvard.edu/abs/2016A%26A...595A...1G/abstract>`_ for the details on OBMT.

    Parameters
    ----------

    release : str
        Data release. One of dr1, dr2, dr3, dr4, dr5.
    epsl : boolean
        Include the Ecliptic Pole Scanning  Law (EPSL) period (default True)
    absolute : boolean
        Return absolute time as astropy.time.Time objects, if false return times as Julian days referred to J2010.0.

    Returns
    -------
    start, end : tuple
        Start and end times in as astropy.time.Time instances, or in Julian days refered to J2010.0

    Notes
    -----
    The end time for DR5 is an estimate.
    """
    _check_release(release)
    if epsl:
        if absolute:
            return _dpac_times["start"]["tcb"], _dpac_times[release]["tcb"]
        else:
            return (
                _dpac_times["start"]["tcb"].jd - _jdref,
                _dpac_times[release]["tcb"].jd - _jdref,
            )
    else:
        if absolute:
            return _dpac_times["end_epsl"]["tcb"], _dpac_times[release]["tcb"]
        else:
            return (
                _dpac_times["end_epsl"]["tcb"].jd - _jdref,
                _dpac_times[release]["tcb"].jd - _jdref,
            )
