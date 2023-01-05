Uncertainty models
==================

Package `pygaia.errors`

This package provides the tools for simulating the Gaia catalogue uncertainties according to the prescriptions on the `Gaia science performance pages <http://www.cosmos.esa.int/web/gaia/science-performance>`_. As of version 3 of PyGaia the uncertainties are simulated for Gaia (E)DR3, and the future Gaia DR4 and DR5 releases. For Gaia DR1/DR2 performance simulations refer to the older versions of Pygaia.

Astrometric uncertainties
-------------------------

Detailed API: :py:mod:`pygaia.errors.astrometric`

Functions for simulating the astrometric uncertainties according to the recipe outlined in the `astrometric section of the Gaia science performance pages <https://www.cosmos.esa.int/web/gaia/science-performance#astrometric%20performance>`_. Uncertainties are provided as a function of :math:`G` for Gaia data release (E)DR3, DR4, or DR5.

Photometric uncertainties
-------------------------

Detailed API: :py:mod:`pygaia.errors.photometric`

Functions for simulating the photometric uncertainties on :math:`G`, :math:`G_\mathrm{BP}`, or :math:`G_\mathrm{RP}`, according to the recipe outlined in the `photometric section of the Gaia science performance pages <https://www.cosmos.esa.int/web/gaia/science-performance#photometric%20performance>`_. Uncertainties are provided as a function of :math:`G`, :math:`G_\mathrm{BP}`, or :math:`G_\mathrm{RP}`, respectively, for Gaia data release (E)DR3, DR4, or DR5.

Radial velocity uncertainties
-----------------------------

Detailed API: :py:mod:`pygaia.errors.spectroscopic`

Functions for simulating the radial velocity uncertainties according to the recipe outlined in the `spectroscopic section of the Gaia science performance pages <https://www.cosmos.esa.int/web/gaia/science-performance#spectroscopic%20performance>`_. Uncertainties are provided as a function of :math:`G_\mathrm{RVS}`, :math:`T_\mathrm{eff}`, :math:`\log g` for Gaia data release (E)DR3, DR4, or DR5.

Utility functions
-----------------

Detailed API: :py:mod:`pygaia.errors.utils`

Utility function used for the uncertainty modelling.
