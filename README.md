# PyGaia

Python toolkit for basic Gaia data simulation and data manipulation.

> **Warning**
> __2022.11.29__ I'm in the process of modifying the code to version 3, intended for Gaia DR3+ performance simulations. The changes to implement the most recent astrometry, photometry, and radial velocity uncertainty models from the [Gaia science performance pages](https://www.cosmos.esa.int/web/gaia/science-performance) are in place. However the code is not yet fully tested and documented, so for now use at your own risk. If you don't want to do this, use the most recent stable version 2.2 from PyPi.

## Description

PyGaia provides python modules for the simulation of Gaia data and their uncertainties, as well modules for the
manipulation of the Gaia catalogue data. In particular transformations between astrometric observables and phase space
variables are provided as well as transformations between sky coordinate systems. Only (very) basic functionality is
provided. Full blown simulations of Gaia data in all their gory detail requires the Java tools developed by the Gaia
Data Processing and Analysis Consortium (DPAC) in particular its Coordination Unit 2 (CU2).

This toolkit is basically an implementation of the performance models for Gaia
which are publicly available at:
[http://www.cosmos.esa.int/web/gaia/science-performance](http://www.cosmos.esa.int/web/gaia/science-performance). In
addition much of the material in chapter 4 of the book [Astrometry for Astrophysics: Methods, Models, and
 Applications (2012, van Altena et al.)](http://www.cambridge.org/9780521519205) is implemented.

* The code in this package is __not intended for accurate astrometry applications__, such as predicting in detail
 astrometric paths of stars on the sky. 
* Epoch transformation _is_ provided, including the transformation of the astrometric covariance matrix to different
 epochs.
 
## Astropy astrometry modules versus pygaia.astrometry

No competition! It is higly recommended to use the [Astropy](https://www.astropy.org/) facilities for handling
 astrometric data, including transformations from (Cartesian) phase space coordinates to astrometric observables and
  vice versa. See the 
[`astropy.coordinates`](https://docs.astropy.org/en/stable/coordinates/index.html) package. Compared to the 
[`pygaia.astrometry`](./pygaia/astrometry) package this gives you more functionality, the use of units, and much better
 maintained code.
 
The only functionality not (yet) provided in Astropy is the transformation of the covariance matrix of the
astrometric observables to to a different epoch. This is implemented in the class 
[`pygaia.astrometry.coordinates.EpochPropagation`](./pygaia/astrometry/coordinates.py). Epoch transformation as such is
implemented in Astropy as the
 [`apply_space_motion()`](https://docs.astropy.org/en/stable/coordinates/apply_space_motion.html) function of
  the [`SkyCoord`](https://docs.astropy.org/en/stable/api/astropy.coordinates.SkyCoord.html#astropy.coordinates.SkyCoord)
   class.

## Documentation

All classes and methods/functions are documented so use the python help() function to find out more.

## Installation

To install from source.

```bash
git clone https://github.com/agabrown/PyGaia.git
cd PyGaia
python -m pip install .
```

## Dependencies

This package is intended for Python3.

The following python packages are required:

* [numpy](https://numpy.org/)
* [scipy](https://scipy.org/)
* [pandas](https://pandas.pydata.org/)
* [UltraJSON](https://github.com/ultrajson/ultrajson)

For the plotting tools:

* [matplotlib](https://matplotlib.org/)
* [Cartopy](https://scitools.org.uk/cartopy/docs/latest/)

## Attribution

Please acknowledge the Gaia Project Scientist Support Team and the Gaia Data Processing and Analysis Consortium 
(DPAC) if you used this code in your research.

## License

Copyright (c) 2012-2022 Anthony Brown, Leiden University, Gaia Data Processing and Analysis Consortium