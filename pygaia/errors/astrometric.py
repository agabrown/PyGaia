__all__ = ['parallax_uncertainty', 'position_uncertainty', 'proper_motion_uncertainty', 'total_position_uncertainty',
           'total_proper_motion_uncertainty']

import numpy as np

from pygaia.errors.utils import calc_z_plx

# Scaling factors for sky averaged position and proper motion uncertainties. The uncertainties are scaled with respect
# to the parallax uncertainty values. Note that the uncertainties are quoted in true arc terms (using phi*) for the
# longitude like component.
_scaling_for_positions = {'dr3': {'Total': 0.75, 'AlphaStar': 0.80, 'Delta': 0.70},
                          'dr4': {'Total': 0.75, 'AlphaStar': 0.80, 'Delta': 0.70},
                          'dr5': {'Total': 0.75, 'AlphaStar': 0.80, 'Delta': 0.70}}
_scaling_for_proper_motions = {'dr3': {'Total': 0.96, 'AlphaStar': 1.03, 'Delta': 0.89},
                               'dr4': {'Total': 0.54, 'AlphaStar': 0.58, 'Delta': 0.50},
                               'dr5': {'Total': 0.27, 'AlphaStar': 0.29, 'Delta': 0.25}}

# Scaling factors for observation time interval with respect to (E)DR3:
# sqrt(33.12/59) for DR4, sqrt(33.12/119) for DR5, for DR4 and DR5 mission lengths of 60 and 120 months (accounting
# for the first month not being used in the astrometric solution).
# Proper motion precision scale as T**(-1.5). The extra factor 1/t is included in the scaling factors above.
#
# The predictions for DR4 and DR5 are based on the (E)DR3 uncertainties, inflated by a 'science margin' of 10 percent
# (factor 1.1).
_science_margin = 1.1
_t_factor = {'dr3': 1.0, 'dr4': 0.749, 'dr5': 0.527}
_default_release = 'dr4'


def parallax_uncertainty(gmag, release=_default_release):
    """
    Calculate the sky averaged parallax uncertainty as a function of G, for a given Gaia data release.

    Parameters
    ----------
    gmag : float or array
        Value(s) of G-band magnitude.

    release : str
        Specify the Gaia data release for which the performance is to be simulated. 'dr3' -> Gaia (E)DR3,
        'dr4' -> Gaia DR4, 'dr5' -> Gaia DR5. Default is 'dr4'.

    Returns
    -------
    parallax_uncertainty: float or array
        The parallax uncertainty in micro-arcseconds.
    """
    z = calc_z_plx(gmag)
    if release == 'dr3':
        return np.sqrt(40 + 800 * z + 30 * z * z) * _t_factor[release] / _science_margin
    else:
        return np.sqrt(40 + 800 * z + 30 * z * z) * _t_factor[release]


def position_uncertainty(gmag, release=_default_release):
    """
    Calculate the sky averaged position uncertainties from G, for a given Gaia data release.

    NOTE! The uncertainties are for sky positions in the ICRS (i.e., right ascension, declination). make sure your
    simulated astrometry is also on the ICRS.

    Parameters
    ----------
    gmag : float or array
        Value(s) of G-band magnitude.

    release : str
        Specify the Gaia data release for which the performance is to be simulated. 'dr3' -> Gaia (E)DR3,
        'dr4' -> Gaia DR4, 'dr5' -> Gaia DR5. Default is 'dr4'.

    Returns
    -------
    ra_cosdelta_uncertainty, delta_uncertainty : float or array
        The uncertainty in alpha* and the uncertainty in delta, in that order, in micro-arcseconds.
    """
    plx_unc = parallax_uncertainty(gmag, release=release)
    return _scaling_for_positions[release]['AlphaStar'] * plx_unc, _scaling_for_positions[release]['Delta'] * plx_unc


def proper_motion_uncertainty(gmag, release=_default_release):
    """
    Calculate the sky averaged proper motion uncertainties from G, for a given Gaia data release.

    NOTE! The uncertainties are for proper motions in the ICRS (i.e., right ascension, declination). make sure
    your simulated astrometry is also on the ICRS.

    Parameters
    ----------
    gmag : float or array
        Value(s) of G-band magnitude.

    release : str
        Specify the Gaia data release for which the performance is to be simulated. 'dr3' -> Gaia (E)DR3,
        'dr4' -> Gaia DR4, 'dr5' -> Gaia DR5. Default is 'dr4'.

    Returns
    -------
    mualpha_cosdelta_uncertainty, mudelta_uncertainty : float or array
        The uncertainty in mu_alpha* and the uncertainty in mu_delta, in that order, in micro-arcseconds/year.
    """
    plx_unc = parallax_uncertainty(gmag, release=release)
    return _scaling_for_proper_motions[release]['AlphaStar'] * plx_unc, \
           _scaling_for_proper_motions[release]['Delta'] * plx_unc


def total_position_uncertainty(gmag, release=_default_release):
    """
    Calculate the sky averaged total position uncertainty as a function of G and for the given Gaia data release.
    This refers to the semi-major axis of the position error ellipse.

    NOTE! The uncertainties are for proper motions in the ICRS (i.e., right ascension, declination). make sure
    your simulated astrometry is also on the ICRS.

    Parameters
    ----------
    gmag : float or array
        Value(s) of G-band magnitude.

    release : str
        Specify the Gaia data release for which the performance is to be simulated. 'dr3' -> Gaia (E)DR3,
        'dr4' -> Gaia DR4, 'dr5' -> Gaia DR5. Default is 'dr4'.

    Returns
    -------
    position_uncertainty
        The semi-major axis of the position error ellipse in micro-arcseconds.
    """
    plx_unc = parallax_uncertainty(gmag, release=release)
    return _scaling_for_positions[release]['Total'] * plx_unc


def total_proper_motion_uncertainty(gmag, release=_default_release):
    """
    Calculate the sky averaged total proper motion uncertainty as a function of G and for the given Gaia data release.
    This refers to the semi-major axis of the proper motion error ellipse.

    NOTE! The uncertainties are for proper motions in the ICRS (i.e., right ascension, declination). make sure
    your simulated astrometry is also on the ICRS.

    Parameters
    ----------
    gmag : float or array
        Value(s) of G-band magnitude.

    release : str
        Specify the Gaia data release for which the performance is to be simulated. 'dr3' -> Gaia (E)DR3,
        'dr4' -> Gaia DR4, 'dr5' -> Gaia DR5. Default is 'dr4'.

    Returns
    -------
    proper_motion_uncertainty
        The semi-major axis of the proper motion error ellipse in micro-arcseconds/year.
    """
    plx_unc = parallax_uncertainty(gmag, release=release)
    return _scaling_for_proper_motions[release]['Total'] * plx_unc
