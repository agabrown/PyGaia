__all__ = ['parallax_uncertainty_sky_avg', 'positionErrorSkyAvg',
           'properMotionErrorSkyAvg', 'parallax_uncertainty', 'positionError', 'properMotionError']

from numpy import sqrt, sin, array, floor, power
from numpy import isscalar
from pygaia.errors.utils import calc_z_plx
from numpy import genfromtxt
from pkg_resources import resource_stream

_table = resource_stream('pygaia', 'data/errorFactorVariationBetaIncludingNumberOfTransits.txt')
_astrometric_error_factors = genfromtxt(_table,
                                        skip_header=4, skip_footer=1,
                                        names=['sinBeta', 'alphaStar', 'delta', 'parallax', 'muAlphaStar', 'muDelta'])
_num_steps_sinbeta = len(_astrometric_error_factors['sinBeta'])

# Scaling factors for sky averaged position and proper motion uncertainties. The uncertainties are scaled with respect
# to the parallax uncertainty values. Note that the uncertainty are quoted in true arc terms (using phi*) for the
# longitude like component.
_scaling_for_positions = {'Total': 0.743, 'AlphaStar': 0.787, 'Delta': 0.699}
_scaling_for_proper_motions = {'Total': 0.526, 'AlphaStar': 0.556, 'Delta': 0.496}

# The upper band of the parallax uncertainties at the bright end
_parallax_error_max_bright = 14.0

# The nominal mission lifetime
_nominal_life_time = 5.0


def uncertainty_scaling_factor(observable, beta):
    """
    Look up the numerical factors to apply to the sky averaged parallax uncertainty in order to obtain uncertainty
    values for a given astrometric parameter, taking the Ecliptic latitude and the number of transits into
    account.

    Parameters
    ----------
    observable : str
        Name of astrometric observable (one of: alphaStar, delta, parallax, muAlphaStar, muDelta)
    beta : float or array
        Values(s) of the Ecliptic latitude.

    Returns
    -------
    uncertainty_scaling_factor : float or array
        Numerical factors to apply to the uncertainties of the given observable.
    """
    if isscalar(beta):
        index = int(floor(abs(sin(beta)) * _num_steps_sinbeta))
        if index == _num_steps_sinbeta:
            return _astrometric_error_factors[observable][_num_steps_sinbeta - 1]
        else:
            return _astrometric_error_factors[observable][index]
    else:
        indices = array(floor(abs(sin(beta)) * _num_steps_sinbeta), dtype=int)
        indices[(indices == _num_steps_sinbeta)] = _num_steps_sinbeta - 1
        return _astrometric_error_factors[observable][indices]


def uncertainty_scaling_mission_length(extension, p):
    """
    Calculate the factor by which to scale uncertainties for a given Gaia mission extension.

    Parameters
    ----------
    extension : float
        The mission extension in years (a negative extension can be used to make crude performance predictions
        part-way through the mission lifetime).
    p : float
        The power by which the uncertainties scale with time (uncertainty ~ t^p, p=-0.5 for parallax and celestial
        position, and -1.5 for proper motion).

    Returns
    -------
    uncertainty_scaling_factor : float
        The factor by which to scale the uncertainties.
    """
    return power((_nominal_life_time + extension) / _nominal_life_time, p)


def parallax_uncertainty_sky_avg(G, vmini, extension=0.0):
    """
    Calculate the sky averaged parallax uncertainty from G and (V-I).

    Parameters
    ----------
    G : float or array
        Value(s) of G-band magnitude.
    vmini : float or array
        Value(s) of (V-I) colour.

    extension : float
        Add this amount of years to the mission lifetime and scale the uncertainties accordingly. Value can be
        negative for shorter mission spans (early data releases).

    Returns
    -------
    parallax_uncertainty: float or array
        The parallax uncertainty in micro-arcseconds.
    """
    factor = uncertainty_scaling_mission_length(extension, -0.5)
    z = calc_z_plx(G)
    return sqrt(-1.631 + 680.766 * z + 32.732 * z * z) * (0.986 + (1.0 - 0.986) * vmini) * factor


def parallax_min_uncertainty(G, vmini, extension=0.0):
    """
    Calculate the minimum parallax  uncertainty at a given G and (V-I). This corresponds to the sky regions with the
    smallest astrometric uncertainties. At the bright end the parallax uncertainty reaches a floor set by calibration
    systematics.

    Parameters
    ----------
    G : float or array
        Value(s) of G-band magnitude.
    vmini : float or array
        Value(s) of (V-I) colour.

    extension : float
        Add this amount of years to the mission lifetime and scale the uncertainties accordingly. Value can be
        negative for shorter mission spans (early data releases).

    Returns
    -------
    minimum_parallax_uncertainty: float or array
        The minimum parallax uncertainty at this G and (V-I) in micro-arcseconds.
    """
    return _astrometric_error_factors["parallax"].min() * parallax_uncertainty_sky_avg(G, vmini, extension=extension)


def parallax_max_uncertainty(G, vmini, extension=0.0):
    """
    Calculate the maximum parallax uncertainty from G and (V-I). This correspond to the sky regions with the
    largest astrometric uncertainty. At the bright end the parallax uncertainty reaches a floor set by calibration
    systematics.

    Parameters
    ----------
    G : float or array
        Value(s) of G-band magnitude.
    vmini : float or array
        Value(s) of (V-I) colour.

    extension : float
        Add this amount of years to the mission lifetime and scale the uncertainties accordingly. Value can be
        negative for shorter mission spans (early data releases).

    Returns
    -------
    maximum_parallax_uncertainty: float or array
        The maximum parallax uncertainty in micro-arcseconds.
    """
    uncertainties = _astrometric_error_factors["parallax"].max() * parallax_uncertainty_sky_avg(G, vmini, extension=extension)
    indices = (uncertainties < _parallax_error_max_bright)
    uncertainties[indices] = _parallax_error_max_bright
    return uncertainties


def parallax_uncertainty(G, vmini, beta, extension=0.0):
    """
    Calculate the parallax uncertainty from G and (V-I) and the Ecliptic latitude beta of the source. The
    parallax uncertainty is calculated by applying numerical factors to the sky average uncertainty. These factors
    depend on beta and the average number of transits across a source for that value of beta.

    The code implements table 6 from the Gaia science performance pages:
    http://www.rssd.esa.int/SA/GAIA/docs/SciencePerformance/table6.htm

    Parameters
    ----------
    G : float or array
        Value(s) of G-band magnitude.
    vmini : float or array
        Value(s) of (V-I) colour.
    beta : float or array
        Value(s) of the Ecliptic latitude.

    extension : float
        Add this amount of years to the mission lifetime and scale the uncertainties accordingly.

    Returns
    -------
    parallax_uncertainty: float or array
        The parallax uncertainty in micro-arcseconds at the given Ecliptic latitudes.
    """
    return parallax_uncertainty_sky_avg(G, vmini, extension=extension) * uncertainty_scaling_factor('parallax', beta)


def positionErrorSkyAvg(G, vmini, extension=0.0):
    """
    Calculate the sky averaged position uncertainties from G and (V-I).

    NOTE! The uncertainties are for sky positions in the ICRS (i.e., right ascension, declination). make sure your
    simulated astrometry is also on the ICRS.

    Parameters
    ----------
    G : float or array
        Value(s) of G-band magnitude.
    vmini : float or array
        Value(s) of (V-I) colour.

    extension : float
        Add this amount of years to the mission lifetime and scale the uncertainties accordingly.

    Returns
    -------
    ra_cosdelta_uncertainty, delta_uncertainty : float or array
        The uncertainty in alpha* and the uncertainty in delta, in that order, in micro-arcseconds.
    """
    parallax_uncertainty = parallax_uncertainty_sky_avg(G, vmini, extension=extension)
    return _scaling_for_positions['AlphaStar'] * parallax_uncertainty, _scaling_for_positions[
        'Delta'] * parallax_uncertainty


def positionMinError(G, vmini, extension=0.0):
    """
    Calculate the minimum position uncertainties from G and (V-I). These correspond to the sky regions with the
    smallest astrometric uncertainties.

    NOTE! The uncertainties are for sky positions in the ICRS (i.e., right ascension, declination). make sure your
    simulated astrometry is also on the ICRS.

    Parameters
    ----------
    G : float or array
        Value(s) of G-band magnitude.
    vmini : float or array
        Value(s) of (V-I) colour.

    extension : float
        Add this amount of years to the mission lifetime and scale the uncertainties accordingly.

    Returns
    -------
    minimum_ra_cosdelta_uncertainty, minimum_delta_uncertainty : float or array
        The minimum uncertainty in alpha* and the uncertainty in delta, in that order, in micro-arcseconds.
    """
    parallax_uncertainty = parallax_uncertainty_sky_avg(G, vmini, extension=extension)
    return _astrometric_error_factors['alphaStar'].min() * parallax_uncertainty, _astrometric_error_factors[
        'delta'].min() * parallax_uncertainty


def positionMaxError(G, vmini, extension=0.0):
    """
    Calculate the maximum position uncertainties from G and (V-I). These correspond to the sky regions with the
    largest astrometric uncertainties.

    NOTE! The uncertainties are for sky positions in the ICRS (i.e., right ascension, declination). make sure your
    simulated astrometry is also on the ICRS.

    Parameters
    ----------
    G : float or array
        Value(s) of G-band magnitude.
    vmini : float or array
        Value(s) of (V-I) colour.

    extension : float
        Add this amount of years to the mission lifetime and scale the uncertainties accordingly.

    Returns
    -------
    maximum_ra_cosdelta_uncertainty, maximum_delta_uncertainty : float or array
        The maximum uncertainty in alpha* and the uncertainty in delta, in that order, in micro-arcseconds.
    """
    parallax_uncertainty = parallax_uncertainty_sky_avg(G, vmini, extension)
    return _astrometric_error_factors['alphaStar'].max() * parallax_uncertainty, _astrometric_error_factors[
        'delta'].max() * parallax_uncertainty


def positionError(G, vmini, beta, extension=0.0):
    """
    Calculate the position uncertainties from G and (V-I) and the Ecliptic latitude beta of the source.

    NOTE! THE uncertainties ARE FOR SKY POSITIONS IN THE ICRS (I.E., RIGHT ASCENSION, DECLINATION). MAKE SURE YOUR
    SIMULATED ASTROMETRY IS ALSO ON THE ICRS.

    Parameters
    ----------
    G : float or array
        Value(s) of G-band magnitude.
    vmini : float or array
        Value(s) of (V-I) colour.
    beta : float or array
        Value(s) of the Ecliptic latitude.

    extension : float
        Add this amount of years to the mission lifetime and scale the uncertainties accordingly.

    Returns
    -------
    ra_cosdelta_uncertainty, delta_uncertainty : float or array
        The uncertainty in alpha* and the uncertainty in delta, in that order, in micro-arcseconds.
    """
    parallax_uncertainty = parallax_uncertainty_sky_avg(G, vmini, extension=extension)
    return uncertainty_scaling_factor('alphaStar', beta) * parallax_uncertainty, uncertainty_scaling_factor('delta',
                                                                                                            beta) * parallax_uncertainty


def properMotionErrorSkyAvg(G, vmini, extension=0.0):
    """
    Calculate the sky averaged proper motion uncertainties from G and (V-I).

    NOTE! The uncertainties are for proper motions in the ICRS (i.e., right ascension, declination). make sure
    your simulated astrometry is also on the ICRS.

    Parameters
    ----------
    G : float or array
        Value(s) of G-band magnitude.
    vmini : float or array
        Value(s) of (V-I) colour.

    extension : float
        Add this amount of years to the mission lifetime and scale the uncertainties accordingly.

    Returns
    -------
    mualpha_cosdelta_uncertainty, mudelta_uncertainty : float or array
        The uncertainty in mu_alpha* and the uncertainty in mu_delta, in that order, in micro-arcseconds/year.
    """
    factor = uncertainty_scaling_mission_length(extension, -1.5)
    parallax_uncertainty = parallax_uncertainty_sky_avg(G, vmini) * factor
    return _scaling_for_proper_motions['AlphaStar'] * parallax_uncertainty, _scaling_for_proper_motions[
        'Delta'] * parallax_uncertainty


def properMotionMinError(G, vmini, extension=0.0):
    """
    Calculate the minimum proper motion uncertainties from G and (V-I). These correspond to the sky regions with
    the smallest astrometric uncertainties.

    NOTE! The uncertainties are for proper motions in the ICRS (i.e., right ascension, declination). make sure
    your simulated astrometry is also on the ICRS.

    Parameters
    ----------
    G : float or array
        Value(s) of G-band magnitude.
    vmini : float or array
        Value(s) of (V-I) colour.

    extension : float
        Add this amount of years to the mission lifetime and scale the uncertainties accordingly.

    Returns
    -------
    minimum_mualpha_cosdelta_uncertainty, minimum_mudelta_uncertainty : float or array
        The minimum uncertainty in mu_alpha* and the uncertainty in mu_delta, in that order, in micro-arcseconds/year.
    """
    factor = uncertainty_scaling_mission_length(extension, -1.5)
    parallax_uncertainty = parallax_uncertainty_sky_avg(G, vmini) * factor
    return _astrometric_error_factors['muAlphaStar'].min() * parallax_uncertainty, _astrometric_error_factors[
        'muDelta'].min() * parallax_uncertainty


def properMotionMaxError(G, vmini, extension=0.0):
    """
    Calculate the maximum proper motion uncertainties from G and (V-I). These correspond to the sky regions with
    the largest astrometric uncertainties.

    NOTE! The uncertainties are for proper motions in the ICRS (i.e., right ascension, declination). make sure
    your simulated astrometry is also on the ICRS.

    Parameters
    ----------
    G : float or array
        Value(s) of G-band magnitude.
    vmini : float or array
        Value(s) of (V-I) colour.

    extension : float
        Add this amount of years to the mission lifetime and scale the uncertainties accordingly.

    Returns
    -------
    maximum_mualpha_cosdelta_uncertainty, maximum_mudelta_uncertainty : float or array
        The maximum uncertainty in mu_alpha* and the uncertainty in mu_delta, in that order, in micro-arcseconds/year.
    """
    factor = uncertainty_scaling_mission_length(extension, -1.5)
    parallax_uncertainty = parallax_uncertainty_sky_avg(G, vmini) * factor
    indices = (parallax_uncertainty < _parallax_error_max_bright)
    parallax_uncertainty[indices] = _parallax_error_max_bright
    return _astrometric_error_factors['muAlphaStar'].max() * parallax_uncertainty, _astrometric_error_factors[
        'muDelta'].max() * parallax_uncertainty


def properMotionError(G, vmini, beta, extension=0.0):
    """
    Calculate the proper motion uncertainties from G and (V-I) and the Ecliptic latitude beta of the source.

    NOTE! The uncertainties are for proper motions in the ICRS (i.e., right ascension, declination). make sure
    your simulated astrometry is also on the ICRS.

    Parameters
    ----------
    G : float or array
        Value(s) of G-band magnitude.
    vmini : float or array
        Value(s) of (V-I) colour.
    beta : float or array
        Value(s) of the Ecliptic latitude.

    extension : float
        Add this amount of years to the mission lifetime and scale the uncertainties accordingly.

    Returns
    -------
    mualpha_cosdelta_uncertainty, mudelta_uncertainty : float or array
        The uncertainty in mu_alpha* and the uncertainty in mu_delta, in that order, in micro-arcseconds/year.
    """
    factor = uncertainty_scaling_mission_length(extension, -1.5)
    parallax_uncertainty = parallax_uncertainty_sky_avg(G, vmini) * factor
    return uncertainty_scaling_factor('muAlphaStar', beta) * parallax_uncertainty, \
           uncertainty_scaling_factor('muDelta', beta) * parallax_uncertainty


def totalProperMotionErrorSkyAvg(G, vmini, extension=0.0):
    """
    Calculate the sky averaged total proper motion uncertainty from G and (V-I). This refers to the uncertainty on the
    length of the proper motion vector.

    NOTE! The uncertainties are for proper motions in the ICRS (i.e., right ascension, declination). make sure
    your simulated astrometry is also on the ICRS.

    Parameters
    ----------
    G : float or array
        Value(s) of G-band magnitude.
    vmini : float or array
        Value(s) of (V-I) colour.

    extension : float
        Add this amount of years to the mission lifetime and scale the uncertainties accordingly.

    Returns
    -------
    proper_motion_uncertainty
        The uncertainty on the total proper motion in micro-arcseconds/year.
    """
    factor = uncertainty_scaling_mission_length(extension, -1.5)
    parallax_uncertainty = parallax_uncertainty_sky_avg(G, vmini) * factor
    return _scaling_for_proper_motions['Total'] * parallax_uncertainty
