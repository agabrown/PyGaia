.. :changelog:

1.2 (December 2016)
+++++++++++++++++++

- Add method to CoordinateTransformation for the transformation of the full (5x5) covariance matrix of
  the astrometric parameters.

- Add keyword to astrometric errors prediction functions that allows to specify an extended mission
  lifetime.

1.1 (September 2016)
++++++++++++++++++++

- Bug fix in coordinate transformation code.

- Updated photometry and radial velocity error predictions.

- End-of-mission photometry errors, including calibration floor, introduced for the broad band fluxes.

- Photometry errors now include 20% margin for CCD-transit and end-of-mission predictions.

- Example plot of photometry errors fixed.

1.0 (November 2015)
+++++++++++++++++++

- Added numerical constants.

- Improved setup.py, made code compatible with Python3

0.9 (September 2015)
++++++++++++++++++++

- Photometric performance predictions updated to post-launch estimates.

0.83 (July 2015)
+++++++++++++++++++

- Minor error in documentation of astrometryToPhaseSpace() method corrected.

0.82 (July 2015)
+++++++++++++++++++

- Error corrected in transformSkyCoordinateErrors() method. Thanks to Teresa Antoja and Taniya Parikh!

0.81 (June 2015)
+++++++++++++++++++

- Forgot to update changelog for version 0.8

0.8 (June 2015)
+++++++++++++++++++

- Radial velocity performance predictions updated to post-launch estimates.

0.7 (December 2014)
+++++++++++++++++++

- Astrometry performance predictions updated to post-launch estimates.

0.6 (July 2014)
+++++++++++++++

- Warning on upcoming changes in performance predictions, following the Gaia
  commissioning period
- radial velocity horizons plot in examples folder

0.5 (August 2013)
+++++++++++++++++

- Utilities for obtaining absolute magnitudes of stars in V and G.
- Functions to obtain the upper and lower bounds on the astrometric parameter
  errors (corresponding to the sky regions with best/worst astrometric errors).
- Proper motion error plot.
- Parallax horizon plot.

0.4 (April 2013)
++++++++++++++++

- Added transformation of proper motions and of position and proper motion errors.

0.31 (February 2013)
++++++++++++++++++++

- Updated README. TODO added.

0.3 (February 2013)
+++++++++++++++++++

- Added documentation on installation requirements. Added the handling of an
  ImportError for the argparse module to the example scripts.

0.2 (February 2013)
+++++++++++++++++++

- Problems in setup.py fixed as well is bugs in the error simulation code.

0.1 (February 2013)
+++++++++++++++++++

- First release

0.0 (October 2012)
++++++++++++++++++

- Creation from bits and pieces of python code that AB had lying around.
