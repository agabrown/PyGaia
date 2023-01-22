"""
Provides transformations between Gaia and a few other photometric systems.

The transformations are taken from the `Gaia DR3 documentation. <https://gea.esac.esa.int/archive/documentation/GDR3/Data_processing/chap_cu5pho/cu5pho_sec_photSystem/cu5pho_ssec_photRelations.html>`_
"""
import numpy as np

__all__ = [
    "gbrminv_from_vminic",
    "gminic_from_vminic",
    "grvs_from_gmingrp",
    "grvs_from_vminic",
]


def gbrminv_from_vminic(vminic):
    r"""
    Transformation from Johnson-Cousins :math:`V-I_\mathrm{c}` to Gaia magnitudes.

    Calculate :math:`G-V`, :math:`G_\mathrm{BP}-V`, :math:`G_\mathrm{RP}-V`, and
    :math:`G_\mathrm{BP}-G_\mathrm{RP}` from the input values of
    :math:`V-I_\mathrm{c}`

    Parameters
    ----------
    vminic : ndarray, float
        Value(s) of :math:`V-I_\mathrm{c}` for which to calculate the Gaia magnitudes

    Returns
    -------
    gminv, gbpminv, grpmin, gbpmingrp : ndarray, float
        The value(s) of :math:`G-V`, :math:`G_\mathrm{BP}-V`,
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


def gminic_from_vminic(vminic):
    r"""
    Transformation from Johnson-Cousins :math:`V-I_\mathrm{c}` to :math:`G-I_\mathrm{c}`.

    Calculate :math:`G-I_\mathrm{c}` from the input values of :math:`V-I_\mathrm{c}`

    Parameters
    ----------
    vminic : ndarray, float
        Value(s) of :math:`V-I_\mathrm{c}` for which to calculate the :math:`G-I_\mathrm{c}` values.

    Returns
    -------
    gminic : ndarray, float
        The value(s) of :math:`G-I_\mathrm{c}`.
    """
    return 0.01753 + 0.76 * vminic - 0.0991 * vminic**2


def grvs_from_gmingrp(gminrp):
    r"""
    Calculate :math:`G_\mathrm{RVS}` from the input value(s) of :math:`G-G_\mathrm{RP}`.
    Use the formulae presented in `Sartoretti et al. (2022)
    <https://ui.adsabs.harvard.edu/abs/2022arXiv220605725S/abstract>`_.

    Parameters
    ----------
    gminrp : ndarray, float
        Value(s) of :math:`G-G_\mathrm{RP}` from which to calculate
        :math:`G_\mathrm{RVS}`. Must be in the range :math:`-0.15\leq
        G-G_\mathrm{RVS}\leq 1.7`

    Returns
    -------
    gminrvs : ndarray, float
        The values of :math:`G_\mathrm{RVS}`. NANs are returned for
        :math:`G-G_\mathrm{RP}` values outside the above quoted range.
    """
    gmrp = np.array(gminrp)
    gminrvs = np.zeros_like(gmrp)
    blue = (gmrp >= -0.15) & (gmrp <= 1.2)
    red = (gmrp > 1.2) & (gmrp <= 1.7)
    gminrvs[blue] = (
        -0.0397
        - 0.2852 * gmrp[blue]
        - 0.0330 * np.power(gmrp[blue], 2)
        - 0.0867 * np.power(gmrp[blue], 3)
    )
    gminrvs[red] = (
        -4.0618
        + 10.0187 * gmrp[red]
        - 9.0532 * np.power(gmrp[red], 2)
        + 2.6089 * np.power(gmrp[red], 3)
    )
    gminrvs[np.logical_not(blue | red)] = np.nan
    return gminrvs


def grvs_from_vminic(vminic):
    r"""
    Calculate :math:`G_\mathrm{RVS}` from the input value(s) of :math:`V-I_\mathrm{c}`.

    Parameters
    ----------
    vminic : ndarray, float
        Value(s) of :math:`V-I_\mathrm{c}` from which to calculate
        :math:`G_\mathrm{RVS}`.

    Returns
    -------
    gminrvs : ndarray, float
        The values of :math:`G_\mathrm{RVS}`. NANs are returned for
        :math:`V-I_\mathrm{c}` values corresponding to :math:`G-G_\mathrm{RP}` values
        outside the valid range listed in `Sartoretti et al. (2022)
        <https://ui.adsabs.harvard.edu/abs/2022arXiv220605725S/abstract>`_.
    """
    gminv, _, grpminv, _ = gbrminv_from_vminic(vminic)
    return grvs_from_gmingrp(gminv - grpminv)
