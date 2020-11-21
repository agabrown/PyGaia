"""
Provides numerical constants useful for astrometric calculations.
"""

import numpy as np

# Astronomical Unit in meter, IAU constant and defining length
au_in_meter = 149597870700.0

# AU expressed in mas*pc or muas*kpc
au_mas_parsec = 1000.0

# Number of seconds in Julian year
julian_year_seconds = 365.25 * 86400.0

# AU expressed in km*yr/s
au_km_year_per_sec = au_in_meter / (julian_year_seconds * 1000.0)

# AU expressed in mas*km*yr/s
au_mas_km_year_per_sec = au_km_year_per_sec * 180 / np.pi * 3600 * 1000

# Definition of parsec
parsec = au_in_meter * 180 * 3600 / np.pi
