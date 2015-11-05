"""
Provides numerical constants useuful for astrometric calculations.
"""

from numpy import pi

# Astronomical Unit in meter, IAU constant and defining length
auInMeter = 149597870700.0

# AU expressed in mas*pc or muas*kpc
auMasParsec = 1000.0

# Number of seconds in Julian year
julianYearSeconds = 365.25 * 86400.0

# AU expressed in km*yr/s
auKmYearPerSec = auInMeter/(julianYearSeconds*1000.0)

# Definition of parsec
parsec = auInMeter * 180*3600/pi
