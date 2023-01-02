r"""
Provides functions for simulating the photometric (:math:`G`, :math:`G_\mathrm{BP}`,
:math:`G_\mathrm{RP}`) uncertainties on the Gaia broad-band photometry for Gaia DR3+.
The code reproduces the uncertainty model described on the `Gaia science performance
<https://www.cosmos.esa.int/web/gaia/science-performance#photometric%20performance>`_
pages.

Code taken (with permission) from the notebook at
https://github.com/gaia-dpci/gaia-dr3-photometric-uncertainties
"""
import os
import pandas as pd
import numpy as np
import scipy.interpolate as interpolate

__all__ = ["LogMagUncertainty", "magnitude_uncertainty"]

_ROOT = os.path.abspath(os.path.dirname(__file__))
_spline_csv_file = os.path.join(_ROOT, "data", "LogErrVsMagSpline2.csv")
_default_release = "dr4"
_nobs_drs_bands = {
    "dr3": {"g": 351, "bp": 40, "rp": 40},
    "dr4": {"g": 620, "bp": 70, "rp": 70},
    "dr5": {"g": 1240, "bp": 140, "rp": 140},
}


class LogMagUncertainty:
    r"""
    Estimate the log(mag) vs mag uncertainty for :math:`G`, :math:`G_\mathrm{BP}`,
    :math:`G_\mathrm{RP}` based on Gaia EDR3 photometry.

    The code in this class is a modified version of the Edr3LogMagUncertainty code in
    the Jypyter notebook EDR3_Photometric_Uncertainties.ipynb at
    https://github.com/gaia-dpci/gaia-dr3-photometric-uncertainties. The additional function estimate_for_maglist() is included to facilitate calculating the uncertainties for the number of observationis corresponding to a given Gaia data release and for the input list of magnitudes.
    """

    def __init__(self):
        """ """
        _df = pd.read_csv(_spline_csv_file)
        splines = dict()
        splines["g"] = self.__init_spline(_df, "knots_G", "coeff_G")
        splines["bp"] = self.__init_spline(_df, "knots_BP", "coeff_BP")
        splines["rp"] = self.__init_spline(_df, "knots_RP", "coeff_RP")
        self.__splines = splines
        self.__nobs = {"g": 200, "bp": 20, "rp": 20}
        self.__nobs_drs = _nobs_drs_bands

    def estimate(
        self, band, nobs: np.array([], int) = 0, mag_range=None, mag_samples=1000
    ):
        """
        Estimate the log(mag) vs mag uncertainty

        Parameters
        ----------
        band : str
            name of the band for which the uncertainties should be estimated (case-insensitive). Must be one of "g", "gbp", or "grp".
        nobs : ndarray, int
            number of observations for which the uncertainties should be estimated.
            Must be a scalar integer value or an array of integer values.
        mag_range : array_like
            Magnitude range over which the spline should be evaluated.
            The default and maximum valid range is (4, 21)
        mag_samples : int
            Number evenly spaced magnitudes (over the mag_range interval) at which the
            splines will be estimated. Default: 1000

        Raises
        ------
        ValueError
            For wrong inputs.

        Returns
        -------
        df : DataFrame
            Pandas dataframe with the interpolated log(mag) uncertainty vs mag.
            The magnitude column is named mag_g, mag_bp, or mag_rp depending of the requested band.
            A column for each value of nobs is provided, in the default case the column is logU_200.
        """
        band = band.lower()
        if band not in ["g", "bp", "rp"]:
            raise ValueError(f"Unknown band: {band}")
        if mag_range is None:
            mag_range = (4.0, 21.0)
        else:
            if mag_range[0] < 4.0:
                raise ValueError(
                    f"Uncertainties can be estimated only in the range {band}[4, 21]"
                )
            elif mag_range[1] > 21.0:
                raise ValueError(
                    f"Uncertainties can be estimated only in the range {band}[4, 21]"
                )
            elif mag_range[0] > mag_range[1]:
                raise ValueError("Malformed magnitude range")
        #
        xx = np.linspace(mag_range[0], mag_range[1], mag_samples)
        __cols = self.__compute_nobs(band, xx, nobs)
        __dc = {f"mag_{band}": xx, **__cols}
        return pd.DataFrame(data=__dc)

    def estimate_for_maglist(
        self, band, maglist: np.array([], float) = 15.0, release=_default_release
    ):
        """
        Estimate the log(mag) vs mag uncertainty

        Parameters
        ----------
        band : str
            name of the band for which the uncertainties should be estimated (case-insensitive). Must be one of "g", "gbp", or "grp".
        maglist : ndarray, float
            List of magnitudes (corresponding to the requested band) for which the
            uncertainties should be estimated. Must be a scalar float value or an array
            of float values. The values must be in the range [4, 21].
        release : str
            Gaia data release for which the uncertainties are simulated (case-insensitive). Must be one of "dr3", "dr4", or "dr5".

        Raises
        ------
        ValueError
            For wrong inputs.

        Returns
        -------
        df : DataFrame
            Pandas dataframe with the interpolated log(mag) uncertainty vs mag.
            The magnitude column is named mag_g, mag_bp, or mag_rp depending of the requested band.
            A column for each value of nobs is provided, in the default case the column is logU_200.
        """
        band = band.lower()
        release = release.lower()
        if band not in ["g", "bp", "rp"]:
            raise ValueError(f"Unknown band: {band}")
        if release not in ["dr3", "dr4", "dr5"]:
            raise ValueError(f"Unknown data release: {release}")
        if np.any(maglist < 4.0) or np.any(maglist > 21.0):
            raise ValueError(
                f"One or more of the values in maglist is outside the range [4, 21]"
            )
        #
        __cols = self.__compute_nobs(band, maglist, self.__nobs_drs[release][band])
        __dc = {f"mag_{band}": maglist, **__cols}
        return pd.DataFrame(data=__dc)

    def __init_spline(self, df, col_knots, col_coeff):
        __ddff = df[[col_knots, col_coeff]].dropna()
        return interpolate.BSpline(
            __ddff[col_knots], __ddff[col_coeff], 3, extrapolate=False
        )

    def __compute_nobs(self, band, xx, nobs):
        if isinstance(nobs, int):
            nobs = [nobs]
        __out = dict()
        for num in nobs:
            if num < 0:
                raise ValueError(f"Number of observations should be strictly positive")
            if num == 0:
                __out[f"logU_{self.__nobs[band]:d}"] = self.__splines[band](xx)
            else:
                __out[f"logU_{num:d}"] = self.__splines[band](xx) - np.log10(
                    np.sqrt(num) / np.sqrt(self.__nobs[band])
                )
        return __out


def magnitude_uncertainty(
    band, maglist: np.array([], float) = 15.0, release=_default_release
):
    r"""
    Provide the uncertainty for :math:`G`, :math:`G_\mathrm{BP}`, and
    :math:`G_\mathrm{RP}` as a function of magnitude (math:`G`, :math:`G_\mathrm{BP}`,
    and :math:`G_\mathrm{RP}`, respectively).

    Parameters
    ----------
    band : str
        name of the band for which the uncertainties are requested (case-insensitive). Must be one of "g", "gbp", or "grp".
    maglist : ndarray, float
        List of magnitudes in the same band for which the uncertainties are requested.
        Must be a scalar float value or an array of float values. The values must be in
        the range [4, 21].
    release : str
        Gaia data release for which the uncertainties are requested (case-insensitive). Must be one of "dr3", "dr4", or "dr5".

    Returns
    -------
    uncs : ndarray, float
        Array with uncertainties in units of mmag.
    """
    lmu = LogMagUncertainty()
    logunc = lmu.estimate_for_maglist(band, maglist, release)[
        f"logU_{_nobs_drs_bands[release.lower()][band.lower()]}"
    ].to_numpy()
    return np.power(10.0, logunc) * 1000
