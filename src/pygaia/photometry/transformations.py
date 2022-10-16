"""
Provides transformations between Gaia and a few other photometric systems.

The transformations are taken from the `Gaia DR3 documentation. <https://gea.esac.esa.int/archive/documentation/GDR3/Data_processing/chap_cu5pho/cu5pho_sec_photSystem/cu5pho_ssec_photRelations.html>`_
"""
import numpy as np

__all__ = ["gbrminv_from_vminic", "gminv_from_vmini", "vmingrvs_from_vmini"]


def gbrminv_from_vminic(vminic):
    r"""
    Transformation from Johnson-Cousins :math:`$V-I_\mathrm{c}$` to Gaia magnitudes.

    Calculate :math:`$G-V$`, :math:`$G_\mathrm{BP}-V$, :math:`G_\mathrm{RP}-V`, and
    :math:`G_\mathrm{BP}-G_\mathrm{RP}` from the input values of
    :math:`$V-I_\mathrm{c}$`

    Parameters
    ----------
    vminic : ndarray, float
        Value(s) of :math:`$V-I_\mathrm{c}$` for which to calculate the Gaia magnitudes

    Returns
    -------
    gminv, gbpminv, grpmin, gbpmingrp : ndarray, float
        The value(s) of :math:`$G-V$`, :math:`$G_\mathrm{BP}-V$,
        :math:`G_\mathrm{RP}-V`, and :math:`G_\mathrm{BP}-G_\mathrm{RP}`
    """
    vminic2 = vminic * vminic
    vminic3 = vminic2 * vminic
    vminic4 = vminic3 * vminic
    gminv = (
        -0.01597
        - 0.02809 * vminic
        - 0.2483 * vminic2
        + 0.03656 * vminic3
        - 0.002939 * vminic4
    )
    gbpminv = -0.0143 + 0.3564 * vminic - 0.1332 * vminic2 + 0.01212 * vminic3
    grpminv = 0.01868 - 0.9028 * vminic - 0.005321 * vminic2 - 0.004186 * vminic3
    gbpmingrp = -0.03298 + 1.259 * vminic - 0.1279 * vminic2 + 0.01631 * vminic3
    return gminv, gbpminv, grpminv, gbpmingrp


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
    return (
        -0.0257 - 0.0924 * vmini - 0.1623 * vmini * vmini + 0.0090 * np.power(vmini, 3)
    )


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
    return (
        0.0119 + 1.2092 * vmini - 0.0188 * vmini * vmini - 0.0005 * np.power(vmini, 3)
    )
