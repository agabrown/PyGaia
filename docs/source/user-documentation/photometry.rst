Photometry utilities
====================

Package: `pygaia.photometry`

This package provides utilities for dealing with Gaia :math:`G`, :math:`G_\mathrm{BP}`, :math:`G_\mathrm{RP}`, and :math:`G_\mathrm{RVS}` photometry. Currently this concerns mainly transformations between Gaia and a few other photometric systems, following the prescriptions in the `Gaia DR3 documentation <https://gea.esac.esa.int/archive/documentation/GDR3/Data_processing/chap_cu5pho/cu5pho_sec_photSystem/cu5pho_ssec_photRelations.html>`_.

Transformations
---------------

Detailed API: :py:mod:`pygaia.photometry.transformations`

rovides transformations between Gaia and a few other photometric systems.

Utilities
---------

Detailed API: :py:mod:`pygaia.photometry.utils`

Provides utility functions for the handling of Gaia photometry. At the moment only functions to calculate photometric colours from the input spectral types.