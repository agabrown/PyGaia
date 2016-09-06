PyGaia
======

**Python toolkit for basic Gaia data simulation, manipulation, and analysis**

PyGaia provides python modules for the simulation of Gaia data and their errors,
as well modules for the manipulation and analysis of the Gaia catalogue data. In
particular transformations between astrometric observables and phase space
variables are provided as well as transformations between sky coordinate
systems. Only (very) basic functionality is provided. Full blown simulations of
Gaia data in all their gory detail requires the Java tools developed by the Gaia
Data Processing and Analysis Consortium (DPAC) in particular its Coordination
Unit 2 (CU2).

This toolkit is basically an implementation of the performance models for Gaia
which are publicly available at:
`<http://www.cosmos.esa.int/web/gaia/science-performance>`_. In
addition much of the material in chapter 4 of the book `Astrometry for
Astrophysics: Methods, Models, and Applications (2012, van Altena et al.)
<http://www.cambridge.org/9780521519205>`_ is implemented.

Note that the code in this package is *not intended for accurate astrometry
applications*, such as predicting in detail astrometric paths of stars on the
sky, or transforming between observation epochs.

Documentation
-------------

All classes and methods/functions are documented so use the python help()
function to find out more. More extensive documentation will follow.

Installation notes
------------------

This package was developed in a python 2.7 environment and you may experience
problems if you have an older version installed. In particular the scripts in
the *examples* folder will not run because they expect the argparse module to be
present.

The following python packages are required:

* `numpy <http://www.numpy.org/>`_
* `scipy <http://www.scipy.org/>`_

For the plotting tools:

* `matplotlib <http://matplotlib.org/>`_
* `Basemap toolkit for matplotlib <http://matplotlib.org/basemap/>`_

Attribution
-----------

Please acknowledge the Gaia Project Scientist Support Team and the Gaia Data
Processing and Analysis Consortium (DPAC) if you used this code in your
research.

License
-------

Copyright (c) 2012-2016 Anthony Brown, Gaia Data Processing and Analysis Consortium

PyGaia is open source and free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your option)
any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program. If not, see `<http://www.gnu.org/licenses/>`_.
