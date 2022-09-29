"""
Provides functions for simulating the radial velocity uncertainties in the Gaia
catalogue data.
"""
import numpy as np

__all__ = ["vrad_error_sky_avg"]

_vradErrorACoeff = {
    "B0V": 0.90,
    "B5V": 0.90,
    "A0V": 1.0,
    "A5V": 1.15,
    "F0V": 1.15,
    "G0V": 1.15,
    "G5V": 1.15,
    "K0V": 1.15,
    "K1IIIMP": 1.15,
    "K4V": 1.15,
    "K1III": 1.15,
}
_vradErrorBCoeff = {
    "B0V": 50.00,
    "B5V": 26.00,
    "A0V": 5.50,
    "A5V": 4.00,
    "F0V": 1.50,
    "G0V": 0.70,
    "G5V": 0.60,
    "K0V": 0.50,
    "K1IIIMP": 0.39,
    "K4V": 0.29,
    "K1III": 0.21,
}
_vradCalibrationFloor = 0.5
_vradMagnitudeZeroPoint = 12.7
_nominal_mission_length = 5.0


def vrad_error_sky_avg(vmag, spt, extension=0.0):
    """
    Calculate radial velocity uncertainty from V and the spectral type. The value of the
    error is an average over the sky.

    Parameters
    ----------
    vmag : float or float array
        Value(s) of V-band magnitude.
    spt : str or str array
        String(s) representing the spectral type of the star.
    extension : float
        Add this amount of years to the mission lifetime and scale the errors
        accordingly. Value can be negative for shorter mission spans (early data
        releases).

    Returns
    -------
    uncertainties : float or float array
        The radial velocity uncertainty in km/s.
    """
    errscaling = 1.0 / np.sqrt(
        (_nominal_mission_length + extension) / _nominal_mission_length
    )
    if np.isscalar(spt):
        return (
            _vradCalibrationFloor
            + _vradErrorBCoeff[spt]
            * np.exp(_vradErrorACoeff[spt] * (vmag - _vradMagnitudeZeroPoint))
            * errscaling
        )
    else:
        uncertainties = np.zeros_like(vmag)
        for s in np.unique(spt):
            s_filt = (spt==s)
            uncertainties[s_filt] = _vradCalibrationFloor + _vradErrorBCoeff[s] * np.exp(
                _vradErrorACoeff[s] * (vmag[s_filt] - _vradMagnitudeZeroPoint)) * errscaling

        return uncertainties
