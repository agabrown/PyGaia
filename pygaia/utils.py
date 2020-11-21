__all__ = ['enum', 'construct_covariance_matrix']

import numpy as np

from pygaia.astrometry.constants import au_km_year_per_sec


def enum(typename, field_names):
    """
    Create a new enumeration type.
  
    Code is copyright (c) Gabriel Genellina, 2010, MIT License.

    Parameters
    ----------

    typename - Name of the enumerated type
    field_names - Names of the fields of the enumerated type
    """

    if isinstance(field_names, str):
        field_names = field_names.replace(',', ' ').split()
    d = dict((reversed(nv) for nv in enumerate(field_names)), __slots__=())
    return type(typename, (object,), d)()


def construct_covariance_matrix(cvec, parallax, radial_velocity, radial_velocity_error):
    """
    Take the astrometric parameter standard uncertainties and the uncertainty correlations as quoted in
    the Gaia catalogue and construct the covariance matrix.

    Parameters
    ----------

    cvec : array_like
        Array of shape (15,) (1 source) or (n,15) (n sources) for the astrometric parameter standard
        uncertainties and their correlations, as listed in the Gaia catalogue [ra_error, dec_error,
        parallax_error, pmra_error, pmdec_error, ra_dec_corr, ra_parallax_corr, ra_pmra_corr,
        ra_pmdec_corr, dec_parallax_corr, dec_pmra_corr, dec_pmdec_corr, parallax_pmra_corr,
        parallax_pmdec_corr, pmra_pmdec_corr]. Units are (mas^2, mas^2/yr, mas^2/yr^2).
    
    parallax : array_like (n elements)
        Source parallax (mas).
    
    radial_velocity : array_like (n elements)
        Source radial velocity (km/s, does not have to be from Gaia RVS!). If the radial velocity is not
        known it can be set to zero.

    radial_velocity_error : array_like (n elements)
        Source radial velocity  uncertainty (km/s). If the radial velocity is not know this can be set to
        the radial velocity dispersion for the population the source was drawn from.

    Returns
    -------

    Covariance matrix as a 6x6 array.
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
        cmat[k, 0:5, 0:5] = cv[k, 0:5] ** 2

    iu = np.triu_indices(5, k=1)
    for k in range(10):
        i = iu[0][k]
        j = iu[1][k]
        cmat[:, i, j] = cv[:, i] * cv[:, j] * cv[:, k + 5]
        cmat[:, j, i] = cmat[:, i, j]

    for k in range(nsources):
        cmat[k, 0:5, 5] = cmat[k, 0:5, 2] * np.atleast_1d(radial_velocity)[k] / au_km_year_per_sec
    cmat[:, 5, 0:5] = cmat[:, 0:5, 5]
    cmat[:, 5, 5] = cmat[:, 2, 2] * (radial_velocity ** 2 + radial_velocity_error ** 2) / au_km_year_per_sec ** 2 + \
                    (parallax * radial_velocity_error / au_km_year_per_sec) ** 2

    return np.squeeze(cmat)
