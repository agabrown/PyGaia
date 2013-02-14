#!/usr/bin/env python

from distribute_setup import use_setuptools
use_setuptools()

try:
  from setuptools import setup
  setup
except ImportError:
  from distutils.core import setup
  setup

import os
import sys

sys.path.append(os.path.dirname(__file__))

import pygaia

extra = {}
if sys.version_info >= (3,):
    extra['use_2to3'] = True

setup(
    name="PyGaia",
    version=pygaia.__version__,
    author="Anthony Brown",
    author_email="brown@strw.leidenuniv.nl",
    url="https://github.com/agabrown/PyGaia",
    license="MIT",
    packages=["pygaia", "pygaia.errors", "pygaia.astrometry", "pygaia.photometry", "pygaia.plot"],
    description="A Python package for basic Gaia data simulation, manipulation, and analysis",
    long_description="\n"+open("README.rst").read() + "\n\n"
                    + "Changelog\n"
                    + "---------\n\n"
                    + open("HISTORY.rst").read(),
    package_data={"": ["LICENSE", "AUTHORS.rst", "HISTORY.rst", "INSTALL"], "pygaia": ["data/*.txt"]},
    include_package_data=True,
    install_requires=[
                        "numpy",
                        "scipy"
                     ],
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU Lesser General Public License v3 or later (LGPLv3+)",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
        "Topic :: Scientific/Engineering :: Astronomy",
    ],
    **extra
)
