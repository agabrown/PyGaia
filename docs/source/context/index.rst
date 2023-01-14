#################
PyGaia in context
#################

Other Gaia catalogue data simulators
====================================

* `Gaia-errors <https://github.com/mromerog/Gaia-errors>`_ Fortran package with functionalities similar to PyGaia.
* `Gaia Object Generator <https://gea.esac.esa.int/archive/documentation/GDR3/Data_processing/chap_cu2sim/sec_cu2GOG/>`_ The DPAC Gaia catalogue simulator. The simulated data is available from the `Gaia Archive <https://gea.esac.esa.int/archive/>`_.
* `Ananke: Gaia on FIRE <https://fire.northwestern.edu/ananke/>`_ Synthetic phase-space surveys designed to resemble Gaia DR2 from different solar viewpoints in the Latte suite of FIRE-2 simulations of Milky Way-mass galaxies

PyGaia and astrometry
=====================

PyGaia is primarily intended as a Python implementation of the `Gaia science performance
predictions <http://www.cosmos.esa.int/web/gaia/science-performance>`_ but does contain
some general tools for working with astrometric data. Such tools can also be found in
Astropy. It is in fact highly recommended to use the `Astropy
<https://www.astropy.org/>`_ facilities for handling astrometric data, including
transformations from (Cartesian) phase space coordinates to astrometric observables and
vice versa. See the `astropy.coordinates
<https://docs.astropy.org/en/stable/coordinates/index.html>`_ package. Compared to the
``pygaia.astrometry`` package this gives you more functionality, the use of units, and
much better, and much better maintained, code.

The only functionality not (yet) provided in Astropy is the propagation of the
covariance matrix of the astrometric observables to to a different epoch. This is
implemented in the class ``pygaia.astrometry.coordinates.EpochPropagation``.  Epoch propagation as such is implemented in Astropy as the `apply_space_motion
<https://docs.astropy.org/en/stable/coordinates/apply_space_motion.html>`_ function of
the `SkyCoord
<https://docs.astropy.org/en/stable/api/astropy.coordinates.SkyCoord.html#astropy.coordinates.SkyCoord>`_
class.

Other Python packages for astrometry
------------------------------------

* `PyERFA <https://github.com/liberfa/pyerfa>`_ Python wrapper for the `ERFA library <https://github.com/liberfa/erfa>`_ (Essential Routines for Fundamental Astronomy), a C library containing key algorithms for astronomy, which is based on the `SOFA library <https://www.iausofa.org/>`_ published by the International Astronomical Union (IAU).
* `pystrometry <https://github.com/Johannes-Sahlmann/pystrometry>`_ Package to support the analysis of high-precision astrometry timeseries, in particular the determination of Keplerian orbits.
* `astromet <https://github.com/zpenoyre/astromet.py>`_  A simple python package for generating astrometric tracks of single stars and the center of light of unresolved binaries, blended and lensed systems. Includes a close emulation of Gaia's astrometric fitting pipeline.
* `HTOF <https://github.com/gmbrandt/HTOF>`_ HTOF: Code which parses the intermediate data from Hipparcos and Gaia and fits astrometric solutions to those data. Capable of computing likelyhoods and parameter errors in line with the catalog.
* `orvara <https://github.com/t-brandt/orvara>`_ Package for fitting orbits of bright stars and their faint companions (exoplanets, brown dwarfs, white dwarfs, and low-mass stars).
* `orbitize <https://orbitize.readthedocs.io/en/latest/>`_ Package for fitting orbits of directly imaged planets.