__all__ = ['calc_z_plx', 'calc_z_gmag', 'calc_z_bprp']

import numpy as np

_bright_floor_star_plx = 13.0
_bright_floor_star_gmag = 12.0
_bright_floor_star_bprp = 11.0


def calc_z_plx(gmag):
    """
    Calculate the value for the parameter z in the formula for parallax errors as a function of G and (V-I).

    Parameters
    ----------
    gmag : float or array
        Value of G-band magnitude.

    Returns
    -------
    z : float or array
        Value of z.
    """
    gatefloor = np.power(10.0, 0.4 * (_bright_floor_star_plx - 15.0))
    if np.isscalar(gmag):
        result = np.amax((gatefloor, np.power(10.0, 0.4 * (gmag - 15.0))))
    else:
        result = np.power(10.0, 0.4 * (gmag - 15.0))
        indices = (result < gatefloor)
        result[indices] = gatefloor
    return result


def calc_z_gmag(gmag):
    """
    Calculate the value for the parameter z in the formula for G magnitude errors as a function of G and (V-I).

    Parameters
    ----------
    gmag : float or array
        Value of G-band magnitude.

    Returns
    -------
    z : float or array
        Value of z.
    """
    gatefloor = np.power(10.0, 0.4 * (_bright_floor_star_gmag - 15.0))
    if np.isscalar(gmag):
        result = np.amax((gatefloor, np.power(10.0, 0.4 * (gmag - 15.0))))
    else:
        result = np.power(10.0, 0.4 * (gmag - 15.0))
        indices = (result < gatefloor)
        result[indices] = gatefloor
    return result


def calc_z_bprp(gmag):
    """
    Calculate the value for the parameter z in the formula for the BP and RP magnitude errors as a
    function of G and (V-I).

    Parameters
    ----------

    gmag : float or array
        Value of G-band magnitude.

    Returns
    -------
    z : float or array
        Value of z for BP/RP.
    """
    gatefloor = np.power(10.0, 0.4 * (_bright_floor_star_bprp - 15.0))
    if np.isscalar(gmag):
        result = np.amax((gatefloor, np.power(10.0, 0.4 * (gmag - 15.0))))
    else:
        result = np.power(10.0, 0.4 * (gmag - 15.0))
        indices = (result < gatefloor)
        result[indices] = gatefloor
    return result
