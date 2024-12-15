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
   :caption: Contents:

   getting-started/index
   context/index
   notebooks/index
   user-documentation/index
   api

Acknowledgements
----------------

PyGaia is based on the effort by Jos de Bruijne to create and maintain the `Gaia Science Performance pages <https://www.cosmos.esa.int/web/gaia/science-performance>`_` (with support from David Katz, Paola Sartoretti, Francesca De Angeli, Dafydd Evans, `Marco Riello <https://github.com/marc0uk>`_, and Josep Manel Carrasco), and benefits from the suggestions and contributions by `Morgan Fouesneau <https://github.com/mfouesneau>`_, `Tom Callingham <https://github.com/TomCallingham>`_, `John Helly <https://github.com/jchelly>`_, `Javier Olivares <https://github.com/olivares-j>`_, `Henry Leung <https://github.com/henrysky>`_, `Johannes Sahlmann <https://github.com/Johannes-Sahlmann>`_.

The photometric uncertainties code in PyGaia is based on the `tool provided by Gaia DPAC <https://www.cosmos.esa.int/web/gaia/dr3-software-tools>`_` to reproduce (E)DR3 Gaia photometric uncertainties described in the `GAIA-C5-TN-UB-JMC-031 <https://dms.cosmos.esa.int/COSMOS/doc_fetch.php?id=1404728>`_` technical note using data presented in `Riello et al (2021) <https://doi.org/10.1051/0004-6361/202039587>`_.

Attribution
-----------

Please acknowledge the Gaia Project Scientist Support Team and the Gaia Data Processing and Analysis Consortium 
(DPAC) if you used this code in your research.

Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
