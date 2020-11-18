__all__ = ['parallaxErrorSkyAvg', 'parallaxErrorSkyAvgAltStartGate', 'positionErrorSkyAvg',
           'properMotionErrorSkyAvg', 'parallaxError', 'positionError', 'properMotionError']

from numpy import sqrt, sin, array, floor, power
from numpy import isscalar
from pygaia.errors.utils import calcZ, calcZAltStartGate
from numpy import genfromtxt
from pkg_resources import resource_stream

_table = resource_stream('pygaia', 'data/errorFactorVariationBetaIncludingNumberOfTransits.txt')
_astrometricErrorFactors = genfromtxt(_table,
                                      skip_header=4, skip_footer=1,
                                      names=['sinBeta', 'alphaStar', 'delta', 'parallax', 'muAlphaStar', 'muDelta'])
_numStepsSinBeta = len(_astrometricErrorFactors['sinBeta'])

# Scaling factors for sky averaged position and proper motion uncertainties. The uncertainties are scaled with respect
# to the parallax uncertainty values. Note that the uncertainty are quoted in true arc terms (using phi*) for the
# longitude like component.
_scalingForPositions = {'Total': 0.743, 'AlphaStar': 0.787, 'Delta': 0.699}
_scalingForProperMotions = {'Total': 0.526, 'AlphaStar': 0.556, 'Delta': 0.496}

# The upper band of the parallax uncertainties at the bright end
_parallaxErrorMaxBright = 14.0

# The nominal mission lifetime
_nominalLifeTime = 5.0


def errorScalingFactor(observable, beta):
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
        index = int(floor(abs(sin(beta)) * _numStepsSinBeta))
        if index == _numStepsSinBeta:
            return _astrometricErrorFactors[observable][_numStepsSinBeta - 1]
        else:
            return _astrometricErrorFactors[observable][index]
    else:
        indices = array(floor(abs(sin(beta)) * _numStepsSinBeta), dtype=int)
        indices[(indices == _numStepsSinBeta)] = _numStepsSinBeta - 1
        return _astrometricErrorFactors[observable][indices]


def errorScalingMissionLength(extension, p):
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
    return power((_nominalLifeTime + extension) / _nominalLifeTime, p)


def parallaxErrorSkyAvg(G, vmini, extension=0.0):
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
    factor = errorScalingMissionLength(extension, -0.5)
    z = calcZ(G)
    return sqrt(-1.631 + 680.766 * z + 32.732 * z * z) * (0.986 + (1.0 - 0.986) * vmini) * factor


def parallaxMinError(G, vmini, extension=0.0):
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
    return _astrometricErrorFactors["parallax"].min() * parallaxErrorSkyAvg(G, vmini, extension=extension)


def parallaxMaxError(G, vmini, extension=0.0):
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
    uncertainties = _astrometricErrorFactors["parallax"].max() * parallaxErrorSkyAvg(G, vmini, extension=extension)
    indices = (uncertainties < _parallaxErrorMaxBright)
    uncertainties[indices] = _parallaxErrorMaxBright
    return uncertainties


def parallaxError(G, vmini, beta, extension=0.0):
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
    return parallaxErrorSkyAvg(G, vmini, extension=extension) * errorScalingFactor('parallax', beta)


def parallaxErrorSkyAvgAltStartGate(G, vmini, extension=0.0):
    """
    Calculate the sky averaged parallax uncertainty from G and (V-I). In this case assume gating starts at G=13.3
    (to simulate bright star worst performance)

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
    parallax_uncertainty : float or array
        The parallax uncertainty in micro-arcseconds.
    """
    factor = errorScalingMissionLength(extension, -0.5)
    z = calcZAltStartGate(G)
    return sqrt(-1.631 + 680.766 * z + 32.732 * z * z) * (0.986 + (1.0 - 0.986) * vmini) * factor


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
    parallax_uncertainty = parallaxErrorSkyAvg(G, vmini, extension=extension)
    return _scalingForPositions['AlphaStar'] * parallax_uncertainty, _scalingForPositions[
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
    parallax_uncertainty = parallaxErrorSkyAvg(G, vmini, extension=extension)
    return _astrometricErrorFactors['alphaStar'].min() * parallax_uncertainty, _astrometricErrorFactors[
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
    parallax_uncertainty = parallaxErrorSkyAvg(G, vmini, extension)
    return _astrometricErrorFactors['alphaStar'].max() * parallax_uncertainty, _astrometricErrorFactors[
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
    parallax_uncertainty = parallaxErrorSkyAvg(G, vmini, extension=extension)
    return errorScalingFactor('alphaStar', beta) * parallax_uncertainty, errorScalingFactor('delta',
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
    factor = errorScalingMissionLength(extension, -1.5)
    parallax_uncertainty = parallaxErrorSkyAvg(G, vmini) * factor
    return _scalingForProperMotions['AlphaStar'] * parallax_uncertainty, _scalingForProperMotions[
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
    factor = errorScalingMissionLength(extension, -1.5)
    parallax_uncertainty = parallaxErrorSkyAvg(G, vmini) * factor
    return _astrometricErrorFactors['muAlphaStar'].min() * parallax_uncertainty, _astrometricErrorFactors[
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
    factor = errorScalingMissionLength(extension, -1.5)
    parallax_uncertainty = parallaxErrorSkyAvg(G, vmini) * factor
    indices = (parallax_uncertainty < _parallaxErrorMaxBright)
    parallax_uncertainty[indices] = _parallaxErrorMaxBright
    return _astrometricErrorFactors['muAlphaStar'].max() * parallax_uncertainty, _astrometricErrorFactors[
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
    factor = errorScalingMissionLength(extension, -1.5)
    parallax_uncertainty = parallaxErrorSkyAvg(G, vmini) * factor
    return errorScalingFactor('muAlphaStar', beta) * parallax_uncertainty, \
           errorScalingFactor('muDelta', beta) * parallax_uncertainty


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
    factor = errorScalingMissionLength(extension, -1.5)
    parallax_uncertainty = parallaxErrorSkyAvg(G, vmini) * factor
    return _scalingForProperMotions['Total'] * parallax_uncertainty
