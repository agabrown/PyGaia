__all__ = ['gminv_from_vmini', 'vmingrvs_from_vmini']

import numpy as np


def gminv_from_vmini(vmini):
    """
    Calculate the value of (G-V) from (V-I).

    Parameters
    ----------
    vmini : float or array
        The value of (V-I).

    Returns
    -------
    gminv : float or array
        The value of (G-V)
    """
    return -0.0257 - 0.0924 * vmini - 0.1623 * vmini * vmini + 0.0090 * np.power(vmini, 3)


def vmingrvs_from_vmini(vmini):
    """
    Calculate (V-Grvs) from (V-I).

    Parameters
    ----------
    vmini : float or array
        The value of (V-I).

    Returns
    -------
    vmingrvs : float or array
        The value of (V-Grvs).
    """
    return 0.0119 + 1.2092 * vmini - 0.0188 * vmini * vmini - 0.0005 * np.power(vmini, 3)
