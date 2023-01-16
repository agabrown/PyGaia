r"""
Provides utility functions for photometric transformations. Values of :math:`M_V` and :math:`V-I_\mathrm{c}` are from `Pickles (1998) <https://ui.adsabs.harvard.edu/abs/1998PASP..110..863P/abstract>`_.
"""
from pygaia.photometry.transformations import gbrminv_from_vminic

__all__ = ["vminic_from_spt", "vabs_from_spt", "gabs_from_spt"]

_spt_to_vminic_vabs = {
    "B0V": (-0.31, -3.5),
    "B1V": (-0.24, -2.7),
    "B5V": (-0.08, 0.0),
    "A0V": (0.01, 0.0),
    "A5V": (0.16, 1.69),
    "F0V": (0.38, 2.98),
    "G0V": (0.67, 4.24),
    "G2V": (0.72, 4.7),
    "G5V": (0.74, 4.78),
    "K0V": (0.87, 5.58),
    "K1IIIMP": (0.99, 1.53),
    "K4V": (1.23, 7.21),
    "K1III": (1.04, 2.16),
    "M0V": (1.71, 8.62),
    "M2V": (2.02, 9.48),
    "M6V": (3.69, 14.2),
    "M0III": (1.65, -0.66),
    "B0I": (-0.22, -7.1),
    "B1I": (-0.16, -6.7),
}


def vminic_from_spt(spt):
    r"""
    Obtain :math:`(V-I_\mathrm{c})` for the input spectral type.

    Parameters
    ----------
    spt : str
        String representing the spectral type of the star.

    Returns
    -------
    vmini : float
        The value of :math:`(V-I_\mathrm{c})`.

    Raises
    ------
    ValueError
        When an unknown spectral type is specified.
    """
    if spt in _spt_to_vminic_vabs:
        return _spt_to_vminic_vabs[spt][0]
    else:
        message = "Unknown spectral type. Allowed values are: "
        for key in _spt_to_vminic_vabs.keys():
            message += key + " "
        raise ValueError(message)


def vabs_from_spt(spt):
    r"""
    Obtain :math:`M_V` (absolute magnitude in :math:`V`-band) for the input spectral type.

    Parameters
    ----------
    spt : str
        String representing the spectral type of the star.

    Returns
    -------
    vabs : float
        The value of :math:`M_V`.

    Raises
    ------
    ValueError
        When an unknown spectral type is specified.
    """
    if spt in _spt_to_vminic_vabs:
        return _spt_to_vminic_vabs[spt][1]
    else:
        message = "Unknown spectral type. Allowed values are: "
        for key in _spt_to_vminic_vabs.keys():
            message += key + " "
        raise ValueError(message)


def gabs_from_spt(spt):
    r"""
    Obtain :math:`M_G` (absolute magnitude in :math:`G`-band) for the input spectral type.

    Parameters
    ----------
    spt : str
        String representing the spectral type of the star.

    Returns
    -------
    gabs : float
        The value of :math:`M_G`.

    Raises
    ------
    ValueError
        When an unknown spectral type is specified.
    """
    if spt in _spt_to_vminic_vabs:
        vminic = vminic_from_spt(spt)
        gminv, _, _, _ = gbrminv_from_vminic(vminic)
        return vabs_from_spt(spt) + gminv
    else:
        message = "Unknown spectral type. Allowed values are: "
        for key in _spt_to_vminic_vabs.keys():
            message += key + " "
        raise ValueError(message)
