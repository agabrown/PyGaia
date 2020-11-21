__all__ = ['g_magnitude_uncertainty', 'bp_magnitude_uncertainty', 'rp_magnitude_uncertainty']

import numpy as np

from pygaia.errors.utils import calc_z_gmag, calc_z_bprp

# Margin to include on predicted standard uncertainties (i.e. multiply prediction by this value).
_science_margin = 1.2

# Mean number of CCDs crossed by a source in the AF field (G-band photometry)
_mean_num_ccds = (7.0 * 9.0 - 1.0) / 7.0

# End-of-mission CCD transit calibration floor on the photometric uncertainties
_eom_calibration_floor_g = 3.0e-3 / np.np.sqrt(_mean_num_ccds)
_eom_calibration_floor_bp = 5.0e-3
_eom_calibration_floor_rp = 5.0e-3


def g_magnitude_uncertainty(gmag):
    """
    Calculate the single-field-of-view-transit photometric standard uncertainty in the G band as a function
    of G. A 20% margin is included.

    Parameters
    ----------
    gmag : float or array
        Value(s) of G-band magnitude.

    Returns
    -------
    gmag_uncertainty : floar or array
        The G band photometric standard uncertainty in units of magnitude.
    """
    z = calc_z_gmag(gmag)
    return 1.0e-3 * np.sqrt(0.04895 * z * z + 1.8633 * z + 0.0001985) * _science_margin


def g_magnitude_uncertainty_eom(gmag, nobs=70, extension=0.0):
    """
    Calculate the end of mission photometric standard uncertainty in the G band as a function
    of G. A 20% margin is included.

    Parameters
    ----------
    gmag : float or array
        Value(s) of G-band magnitude.

    nobs : int
        Number of observations collected (default 70).
    extension : float
        Add this amount of years to the mission lifetime and scale the uncertainties accordingly. Value can be
        negative for shorter mission spans (early data releases).

    Returns
    -------
    gmag_uncertainty : floar or array
        The G band photometric standard uncertainty in units of magnitude.
    """
    nobs_scaled = round((5.0 + extension) / 5.0 * nobs)
    return np.sqrt((np.power(g_magnitude_uncertainty(gmag) / _science_margin, 2) +
                    _eom_calibration_floor_g * _eom_calibration_floor_g) / nobs_scaled) * _science_margin


def bp_magnitude_uncertainty(gmag, vmini):
    """
    Calculate the single-field-of-view-transit photometric standard uncertainty in the BP band as a function
    of G and (V-I). Note: this refers to the integrated flux from the BP spectrophotometer. A margin of 20%
    is included.

    Parameters
    ----------
    gmag : float or array
        Value(s) of G-band magnitude.
    vmini : float or array
        Value(s) of (V-I) colour.

    Returns
    -------
    bpmag_uncertainty : float or array
        The BP band photometric standard uncertainty in units of magnitude.
    """
    z = calc_z_bprp(gmag)
    a = -0.000562 * np.power(vmini, 3) + 0.044390 * vmini * vmini + 0.355123 * vmini + 1.043270
    b = -0.000400 * np.power(vmini, 3) + 0.018878 * vmini * vmini + 0.195768 * vmini + 1.465592
    c = +0.000262 * np.power(vmini, 3) + 0.060769 * vmini * vmini - 0.205807 * vmini - 1.866968
    return 1.0e-3 * np.sqrt(np.power(10.0, a) * z * z + np.power(10.0, b) * z + np.power(10.0, c))


def bp_magnitude_uncertainty_eom(gmag, vmini, nobs=70, extension=0.0):
    """
    Calculate the end-of-mission photometric standard uncertainty in the BP band as a function of G and (V-I).
    Note: this refers to the integrated flux from the BP spectrophotometer. A margin of 20% is included.

    Parameters
    ----------

    gmag : float or array
        Value(s) of G-band magnitude.
    vmini : float or array
        Value(s) of (V-I) colour.

    nobs : int
        Number of observations collected (default 70).
    extension : float
        Add this amount of years to the mission lifetime and scale the uncertainties accordingly. Value can be
        negative for shorter mission spans (early data releases).

    Returns
    -------
    bpmag_uncertainty : float or array
        The BP band photometric standard uncertainty in units of magnitude.
    """
    nobs_scaled = round((5.0 + extension) / 5.0 * nobs)
    return np.sqrt((np.power(bp_magnitude_uncertainty(gmag, vmini) / _science_margin, 2) +
                    _eom_calibration_floor_bp * _eom_calibration_floor_bp) / nobs_scaled) * _science_margin


def rp_magnitude_uncertainty(gmag, vmini):
    """
    Calculate the single-field-of-view-transit photometric standard uncertainty in the RP band as a function
    of G and (V-I). Note: this refers to the integrated flux from the RP spectrophotometer. A margin of 20%
    is included.

    Parameters
    ----------
    gmag : float or array
        Value(s) of G-band magnitude.
    vmini : float or array
        Value(s) of (V-I) colour.

    Returns
    -------
    rpmag_uncertainty : float or array
        The RP band photometric standard uncertainty in units of magnitude.
    """
    z = calc_z_bprp(gmag)
    a = -0.007597 * np.power(vmini, 3) + 0.114126 * vmini * vmini - 0.636628 * vmini + 1.615927
    b = -0.003803 * np.power(vmini, 3) + 0.057112 * vmini * vmini - 0.318499 * vmini + 1.783906
    c = -0.001923 * np.power(vmini, 3) + 0.027352 * vmini * vmini - 0.091569 * vmini - 3.042268
    return 1.0e-3 * np.sqrt(np.power(10.0, a) * z * z + np.power(10.0, b) * z + np.power(10.0, c))


def rp_magnitude_uncertainty_eom(gmag, vmini, nobs=70, extension=0.0):
    """
    Calculate the end-of-mission photometric standard uncertainty in the RP band as a function of G and (V-I).
    Note: this refers to the integrated flux from the RP spectrophotometer. A margin of 20% is included.

    Parameters
    ----------
    gmag : float or array
        Value(s) of G-band magnitude.
    vmini : float or array
        Value(s) of (V-I) colour.

    nobs : int
        Number of observations collected (default 70).
    extension : float
        Add this amount of years to the mission lifetime and scale the uncertainties accordingly. Value can be
        negative for shorter mission spans (early data releases).

    Returns
    -------
    rpmag_uncertainty : float or array
        The RP band photometric standard uncertainty in units of magnitude.
    """
    nobs_scaled = round((5.0 + extension) / 5.0 * nobs)
    return np.sqrt((np.power(rp_magnitude_uncertainty(gmag, vmini) / _science_margin, 2) +
                    _eom_calibration_floor_rp * _eom_calibration_floor_rp) / nobs_scaled) * _science_margin
