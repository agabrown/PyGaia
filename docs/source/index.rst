.. PyGaia documentation master file, created by
   sphinx-quickstart on Sat Sep 24 18:55:47 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

PyGaia
======

Python modules for the simulation and basic manipulation of Gaia catalogue data and the
corresponding uncertainties. To simulate Gaia astrometric data the following functionalities are provided:

* transform phase space variables to astrometric observables as well as radial velocities and vice versa
* transformations between sky coordinate systems
* epoch transformations of the astrometric data, including the transformation of the astrometric covariance matrix
* uncertainty models for the astrometric, photometric, and radial velocity data from Gaia

This toolkit is basically an implementation of the performance models for Gaia which are
publicly available at: `http://www.cosmos.esa.int/web/gaia/science-performance
<http://www.cosmos.esa.int/web/gaia/science-performance>`_. In addition much of the
material in chapter 4 of the book `Astrometry for Astrophysics: Methods, Models, and
Applications (2012, van Altena et al.) <http://www.cambridge.org/9780521519205>`_ is
implemented.

.. warning:: 
   The code in this package is intended for simulating Gaia catalogue data and its uncertainties and manipulating (Gaia) astrometric data, but **is not intended for accurate on-sky astrometry applications**, such as predicting in detail astrometric paths of stars on the sky.
   
.. toctree::
   :maxdepth: 1

   getting-started/index
   context/index
   user-documentation/index
   api

Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
