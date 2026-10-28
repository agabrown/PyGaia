# PyGaia

[![PyPI](https://img.shields.io/pypi/v/pygaia.svg)](https://pypi.org/project/PyGaia/)
[![Python versions](https://img.shields.io/pypi/pyversions/pygaia.svg)](https://pypi.org/project/PyGaia/)
[![CI](https://github.com/agabrown/PyGaia/actions/workflows/ci.yaml/badge.svg?branch=master)](https://github.com/agabrown/PyGaia/actions/workflows/ci.yaml)
[![Documentation](https://readthedocs.org/projects/pygaia/badge/?version=latest)](https://pygaia.readthedocs.io/en/latest/)
[![License: BSD-3-Clause](https://img.shields.io/badge/license-BSD--3--Clause-blue.svg)](LICENSE)

Python toolkit for Gaia science performance simulation and astrometric catalogue data manipulation.

## Description

PyGaia provides python modules for the simulation of Gaia data and their uncertainties, as well modules for the
manipulation of the Gaia catalogue data. In particular transformations between astrometric observables and phase space
variables are provided as well as transformations between sky coordinate systems. Only (very) basic functionality is
provided. Full blown simulations of Gaia data in all their gory detail requires the Java tools developed by the Gaia
Data Processing and Analysis Consortium (DPAC) in particular its Coordination Unit 2 (CU2).

This toolkit is basically an implementation of the performance models for Gaia which are publicly available at:
[http://www.cosmos.esa.int/web/gaia/science-performance](http://www.cosmos.esa.int/web/gaia/science-performance). In
addition much of the material in chapter 4 of the book [Astrometry for Astrophysics: Methods, Models, and
 Applications (2012, van Altena et al.)](http://www.cambridge.org/9780521519205) is implemented.

* The code in this package is __not intended for accurate astrometry applications__, such as predicting in detail
 astrometric paths of stars on the sky. 
* Epoch transformation _is_ provided, including the transformation of the astrometric covariance matrix to different
 epochs.
 
## Astropy astrometry modules versus pygaia.astrometry

It is recommended to use the [Astropy](https://www.astropy.org/) facilities for handling astrometric data, including transformations from (Cartesian) phase space coordinates to astrometric observables and vice versa. See the [`astropy.coordinates`](https://docs.astropy.org/en/stable/coordinates/index.html) package. Compared to the [`pygaia.astrometry`](./src/pygaia/astrometry) package this gives you more functionality, the use of units, and much better maintained code.
 
The only functionality not (yet) provided in Astropy is the propagation of the covariance matrix of the astrometric observables to to a different epoch. This is implemented in the class [`pygaia.astrometry.coordinates.EpochPropagation`](./src/pygaia/astrometry/coordinates.py). Epoch propagation as such is implemented in Astropy as the [`apply_space_motion()`](https://docs.astropy.org/en/stable/coordinates/apply_space_motion.html) function of the [`SkyCoord`](https://docs.astropy.org/en/stable/api/astropy.coordinates.SkyCoord.html#astropy.coordinates.SkyCoord) class.

## Documentation

The documentation, including example notebooks, is available at [pygaia.readthedocs.io](https://pygaia.readthedocs.io/en/latest/). All classes and methods/functions are also documented, so you can use the python `help()` function to find out more.

## Installation

To install from source.

```bash
git clone https://github.com/agabrown/PyGaia.git
cd PyGaia
python -m pip install .
```

From PyPI:

```bash
pip install pygaia
```

For development, install the package in editable mode with the test dependencies:

```bash
python -m pip install -e ".[test]"
```

See [CONTRIBUTING.md](CONTRIBUTING.md) for more details.

## Dependencies

PyGaia requires Python 3.10 or later.

The following python packages are required:

* [numpy](https://numpy.org/)
* [scipy](https://scipy.org/)
* [pandas](https://pandas.pydata.org/)
* [Astropy](https://www.astropy.org/)

For the plotting tools:

* [matplotlib](https://matplotlib.org/)

## Acknowledgements

PyGaia is based on the effort by Jos de Bruijne to create and maintain the [Gaia Science Performance pages](https://www.cosmos.esa.int/web/gaia/science-performance) (with support from David Katz, Paola Sartoretti, Francesca De Angeli, Dafydd Evans, [Marco Riello](https://github.com/marc0uk), and Josep Manel Carrasco), and benefits from the suggestions and contributions by [Morgan Fouesneau](https://github.com/mfouesneau), [Tom Callingham](https://github.com/TomCallingham), [John Helly](https://github.com/jchelly), [Javier Olivares](https://github.com/olivares-j), [Henry Leung](https://github.com/henrysky), [Johannes Sahlmann](https://github.com/Johannes-Sahlmann).

The photometric uncertainties code in PyGaia is based on the [tool provided by Gaia DPAC](https://www.cosmos.esa.int/web/gaia/dr3-software-tools) to reproduce (E)DR3 Gaia photometric uncertainties described in the [GAIA-C5-TN-UB-JMC-031](https://dms.cosmos.esa.int/COSMOS/doc_fetch.php?id=1404728) technical note using data presented in [Riello et al (2021)](https://doi.org/10.1051/0004-6361/202039587).

## Contributors

Direct contributions to the code base by [John Helly](https://github.com/jchelly), [Javier Olivares](https://github.com/olivares-j), [Tom Callingham](https://github.com/TomCallingham), [Henry Leung](https://github.com/henrysky), [Morgan Fouesneau](https://github.com/mfouesneau), and [Johannes Sahlmann](https://github.com/Johannes-Sahlmann). Contributions are welcome, see [CONTRIBUTING.md](CONTRIBUTING.md).

## Citation

If you use PyGaia in your research, please cite it using the metadata in [CITATION.cff](CITATION.cff) (GitHub shows a "Cite this repository" button in the sidebar). Please also acknowledge the Gaia Project Scientist Support Team and the Gaia Data Processing and Analysis Consortium (DPAC).

## License

Copyright (c) 2012-2024 Anthony Brown, Leiden University, Gaia Data Processing and Analysis Consortium