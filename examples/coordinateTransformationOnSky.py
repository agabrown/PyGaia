"""
Example usage of the plotCoordinateTransformationOnSky() function which shows
how coordinate systems such as ICRS and Galactic are related on the sky.

Anthony Brown May 2019
"""

from pygaia.astrometry.coordinates import Transformations
from pygaia.plot.sky import plot_coordinate_transformation_on_sky

# Ecliptic coordinates on Galactic coordinate sky
plot_coordinate_transformation_on_sky(Transformations.ECL2GAL)
