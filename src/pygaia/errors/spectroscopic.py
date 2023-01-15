"""
Provides functions for simulating the radial velocity uncertainties in the Gaia DR3+
catalogue data. The code reproduces the uncertainty model described on the Gaia `science
performance pages
<https://www.cosmos.esa.int/web/gaia/science-performance#spectroscopic%20performance>`_.
"""
import numpy as np
import ujson
import os

__all__ = ["radial_velocity_uncertainty"]

_default_release = "dr4"

_ROOT = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(_ROOT, "data", "rv_uncertainty_model_coeffs.json")) as fp:
    _rv_unc_model = ujson.load(fp)
    fp.close()

_rv_nb_transits_dr4 = 32
_rv_nb_transits_dr5 = 64
_grvs_zeropoint = 21.317
_exposure_time = 4.4167032
_n_spectrum_per_transit = 3
_pixel_width_al = 0.02453
_band_width = 24.0
_median_background = 4.7
_n_ac_samples_bright = 10
_n_ac_samples_faint = 1
_ron = 3.2


def _in_interval(a: np.array([], float), left, right, closed="both"):
    """
    Check whether the input number or array elements lie within the given interval.

    Parameters
    ----------
    a : ndarray, float
        Input number(s) to check.
    left : float
        Left bound of interval, where left<=right.
    right : float
        Right bound of interval, where right>=left.
    closed : str
        Can be 'both' for a closed interval, 'neither' for an open interval, 'left' for a left-closed interval, or 'right' for a right-closed interval.

    Returns
    -------
    result : ndarray, boolean
        Result of the check, True or False
    """
    if left > right:
        raise ValueError("Left must be less than or equal to right")
    if closed == "both":
        return (a >= left) & (a <= right)
    elif closed == "left":
        return (a >= left) & (a < right)
    elif closed == "right":
        return (a > left) & (a <= right)
    elif closed == "neither":
        return (a > left) & (a < right)
    else:
        raise ValueError("The closed parameter must be one of both|left|right|neither")


def radial_velocity_uncertainty(
    grvs: np.array([], float),
    teff: np.array([], float),
    logg: np.array([], float),
    release=_default_release,
):
    r"""
    Simulate the Gaia DR3 radial velocity uncertainty for the input list of
    :math:`G_\mathrm{RVS}`, :math:`T_\mathrm{eff}`, and :math:`\log(g)` values.

    Parameters
    ----------
    grvs : ndarray, float
        Value(s) of :math:`G_\mathrm{RVS}` for which the calculate the radial velocity uncertainty.
    teff : ndarray, float
        Value(s) of :math:`T_\mathrm{eff}` (in K) for which to calculate the radial velocity uncertainty.
    logg : ndarray, float
        Value(s) of :math:`\log(g)` for which to calculate the radial velocity uncertainty.

    Returns
    -------
    sigma_rv : ndarray, float
        Value(s) of the radial velocity uncertainty. NaNs are returned for input
        magnitudes, temperatures, and surface gravities outside the grids as defined in
        `Katz et al (2022)
        <https://ui.adsabs.harvard.edu/abs/2022arXiv220605902K/abstracti>`_ (their
        figures E.1 and F.1), or outside the model validity ranges quoted on the `Gaia
        Science Performance pages
        <https://www.cosmos.esa.int/web/gaia/science-performance#spectroscopic%20performance>`_.
    """
    ggrvs = np.array(grvs)
    tteff = np.array(teff)
    llogg = np.array(logg)
    if not (
        (tteff.size == llogg.size)
        and (tteff.size == ggrvs.size)
        and (llogg.size == ggrvs.size)
    ):
        raise ValueError("Arrays grvs, teff. and logg must be of same size")
    rv_uncs = np.full(tteff.shape, np.nan)
    if release.upper() == "DR3":
        coeffs = "dr3coeffs"
    else:
        coeffs = "dr45coeffs"
        if release.upper() == "DR4":
            rv_nb_transits = _rv_nb_transits_dr4
        else:
            rv_nb_transits = _rv_nb_transits_dr5

    for group in _rv_unc_model.keys():
        slots = _in_interval(
            tteff,
            _rv_unc_model[group]["teff"][0],
            _rv_unc_model[group]["teff"][1],
            _rv_unc_model[group]["teff"][2],
        ) & _in_interval(
            llogg,
            _rv_unc_model[group]["logg"][0],
            _rv_unc_model[group]["logg"][1],
            _rv_unc_model[group]["logg"][2],
        )
        if np.any(slots):
            if release.upper() == "DR3":
                rv_uncs[slots] = _rv_unc_model[group][coeffs]["sfloor"] + _rv_unc_model[
                    group
                ][coeffs]["b"] * np.exp(
                    _rv_unc_model[group][coeffs]["a"]
                    * (ggrvs[slots] - _rv_unc_model[group][coeffs]["grvs0"])
                )
            else:
                rv_expected_sig_to_noise = np.zeros_like(ggrvs)
                collected_signal = (
                    np.power(10.0, 0.4 * (_grvs_zeropoint - ggrvs))
                    * _exposure_time
                    * _n_spectrum_per_transit
                    * rv_nb_transits
                    * (_pixel_width_al / _band_width)
                )
                bck_per_sample = (
                    _median_background * _n_spectrum_per_transit * rv_nb_transits
                )
                rn_per_sample = _ron * _ron * _n_spectrum_per_transit * rv_nb_transits
                bright = ggrvs <= 7
                faint = np.logical_not(bright)
                rv_expected_sig_to_noise[bright] = collected_signal[bright] / np.sqrt(
                    collected_signal[bright]
                    + bck_per_sample * _n_ac_samples_bright
                    + rn_per_sample * _n_ac_samples_bright
                )
                rv_expected_sig_to_noise[faint] = collected_signal[faint] / np.sqrt(
                    collected_signal[faint]
                    + bck_per_sample * _n_ac_samples_faint
                    + rn_per_sample * _n_ac_samples_faint
                )
                slowsnr = _rv_unc_model[group][coeffs]["sbreak"] * np.power(
                    rv_expected_sig_to_noise[slots]
                    / _rv_unc_model[group][coeffs]["snrbreak"],
                    _rv_unc_model[group][coeffs]["f"],
                )
                shighsnr = _rv_unc_model[group][coeffs]["sfloor"] + (
                    _rv_unc_model[group][coeffs]["sbreak"]
                    - _rv_unc_model[group][coeffs]["sfloor"]
                ) * np.exp(
                    _rv_unc_model[group][coeffs]["g"]
                    * (
                        np.log10(rv_expected_sig_to_noise[slots])
                        - np.log10(_rv_unc_model[group][coeffs]["snrbreak"])
                    )
                )
                h = (
                    1
                    + np.tanh(
                        _rv_unc_model[group][coeffs]["k"]
                        * (
                            np.log10(rv_expected_sig_to_noise[slots])
                            - np.log10(_rv_unc_model[group][coeffs]["snrbreak"])
                        )
                    )
                ) / 2
                rv_uncs[slots] = h * shighsnr + (1 - h) * slowsnr

    if release.upper() == "DR3":
        rv_uncs[(ggrvs > 14) | (rv_uncs > 20.0)] = np.nan
    else:
        rv_uncs[(ggrvs > 16)] = np.nan
        rv_uncs[(ggrvs > 12) & (tteff > 7000)] = np.nan

    return rv_uncs
