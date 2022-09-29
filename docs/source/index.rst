.. PyGaia documentation master file, created by
   sphinx-quickstart on Sat Sep 24 18:55:47 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

PyGaia
======

Python modules for the simulation and basic manipulation of Gaia catalogue data and the
corresponding uncertainties. In particular the following transformations are provided:

* astrometric observables to phase space variables and vice versa
* transformations between sky coordinate systems
* epoch transformations of the astrometric data, including the transformation of the astrometric covariance matrix

This toolkit is basically an implementation of the performance models for Gaia which are
publicly available at: `http://www.cosmos.esa.int/web/gaia/science-performance
<http://www.cosmos.esa.int/web/gaia/science-performance>`_.  In addition much of the
material in chapter 4 of the book `Astrometry for Astrophysics: Methods, Models, and
Applications (2012, van Altena et al.) <http://www.cambridge.org/9780521519205>`_ is
implemented.

.. warning:: 
   The code in this package is **not intended for accurate astrometry applications**,
   such as predicting in detail astrometric paths of stars on the sky.
   
.. toctree::
   :maxdepth: 1

   getting-started/index
   context/index
   api

Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
