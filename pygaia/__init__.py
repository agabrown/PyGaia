"""
Basic Gaia data simulation, manipulation, and analysis toolkit.

PyGaia provides python modules for the simulation of Gaia data and their errors, as well modules for the
manipulation and analysis of the Gaia catalogue data. In particular transformations between astrometric
observables and phase space variables are provided as well as transformations between sky coordinate
systems. Only (very) basic functionality is provided. Full blown simulations of Gaia data in all their
gory detail requires the Java tools developed by the Gaia Data Processing and Analysis Consortium (DPAC)
in particular its Coordination Unit 2 (CU2).

pygaia requires numpy and scipy.
"""

__version__ = "1.2"

try:
    import numpy
except ImportError:
    raise ImportError('NumPy does not seem to be installed.')

try:
    import scipy
except ImportError:
    raise ImportError('SciPy does not seem to be installed.')

__modules__ = ['utils']
