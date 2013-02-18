import os
import sys
import re

from distribute_setup import use_setuptools
use_setuptools()

try:
  from setuptools import setup
  setup
except ImportError:
  from distutils.core import setup
  setup

if sys.version_info >= (3,):
    extra['use_2to3'] = True

vre = re.compile("__version__ = \"(.*?)\"")
m = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "pygaia", "__init__.py")).read()
pygaiaVersion = vre.findall(m)[0]

setup(
    name="PyGaia",
    description="Basic Gaia data simulation, manipulation, and analysis toolkit",
    version=pygaiaVersion,
    author="Anthony Brown",
    author_email="brown@strw.leidenuniv.nl",
    url="https://github.com/agabrown/PyGaia",
    license="LGPLv3+",
    long_description="\n"+open("README.rst").read() + "\n\n"
    + "Changelog\n"
    + "---------\n\n"
    + open("HISTORY.rst").read(),
    packages=['pygaia', 'pygaia.errors', 'pygaia.astrometry', 'pygaia.photometry', 'pygaia.plot',
      'pygaia.tests'],
    package_data={'': ['LICENSE', 'AUTHORS.rst', 'HISTORY.rst', 'INSTALL', 'MANIFEST.in'],
      'pygaia': ['data/*.txt']},
    include_package_data=True,
    install_requires=[
      "numpy",
      "scipy"
      ],
    classifiers=[
      "Development Status :: 3 - Alpha",
      "Intended Audience :: Developers",
      "Intended Audience :: Science/Research",
      "License :: OSI Approved :: GNU Lesser General Public License v3 or later (LGPLv3+)",
      "Operating System :: OS Independent",
      "Programming Language :: Python",
      "Topic :: Scientific/Engineering :: Astronomy",
      ],
    )
