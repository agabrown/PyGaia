Astrometry tools
================

Package: `pygaia.astrometry`

This package provides functions for transforming between astrometric and (Cartesian)
phase space data, the transformations between different coordinate systems, and the epoch propagation of astrometric data. Vector astrometry tools are also included.

Constants
---------

Detailed API: :py:mod:`pygaia.astrometry.constants`

This module contains a number of useful constants for handling astrometric data.

Coordinates
-----------

Detailed API: :py:mod:`pygaia.astrometry.coordinates`

This module provides sky coordinate transformations, epoch propagation, and a function to calculate the angular distance between points on the sky.

* Coordinate transformations between ICRS, Ecliptic, and Galactic coordinates (either sky or phase space coordinates). Note that this can be achieved also with the `astropy coordinates <https://docs.astropy.org/en/stable/coordinates/index.html>`_ package.
* Transformations of uncertainties and covariance matrices from one coordinate system to another are also provided. 
* Epoch propagation from a reference to a future or past epoch, accounting for proper motion and radial velocity. This includes the propagation of the covariance matrix of the astrometric and radial velocity observables. The Astropy `apply_space_motion <https://docs.astropy.org/en/stable/coordinates/apply_space_motion.html>`_ function of the `SkyCoord <https://docs.astropy.org/en/stable/api/astropy.coordinates.SkyCoord.html#astropy.coordinates.SkyCoord>`_ class provides the same functionality but does not include the covariance matrix propagation.

An example coordinate transformation (from Galactic to Ecliptic):

.. code-block:: python
    
    import numpy as np
    from pygaia.astrometry.coordinates import CoordinateTransformation, Transformations
    ct = CoordinateTransformation(Transformations.GAL2ECL)
    
    l = np.deg2rad(30.0)
    b = np.deg2rad(-70.0)
    pml = 5.0
    pmb = -3.0

    lam, beta = ct.transform_sky_coordinates(l, b)
    pmlam, pmbeta = ct.transform_proper_motions(l, b, pml, pmb)

The sky coordinate transformations can be visualized with the functions in the :py:mod:`pygaia.plot.sky` module.

Vector astrometry
-----------------

Detailed API: :py:mod:`pygaia.astrometry.vectorastrometry`

This module provides various vector astrometry functions. These are derived from the material presented in  in chapter 4 of the book `Astrometry for Astrophysics: Methods, Models, and Applications (2012, van Altena et al.) <http://www.cambridge.org/9780521519205>`_. The functionalities provided are:

* Transformations between spherical and Cartesian coordinates.
* Calculation of the `normal triad <https://agabrown.github.io/icrs-coordinates/>`_ from the input sky coordinates.
* Transformations from phase space variables (:math:`x, y, z, v_x, v_y, v_z`) to astrometry and radial velocity (:math:`\alpha, \delta, \varpi, \mu_{\alpha*}, \mu_\delta, v_\mathrm{rad}`)
* Elementary rotation matrices. This is not needed with `SciPy's rotation module <https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.transform.Rotation.html#scipy.spatial.transform.Rotation>`_ available, but I was too lazy to change the code.