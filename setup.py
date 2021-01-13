"""
Setup module for PyGaia.

Anthony Brown 2012--2021

Based on:
https://packaging.python.org/en/latest/distributing.html
https://github.com/pypa/sampleproject
"""

import re
from codecs import open
from os import path

from setuptools import setup, find_packages

this_folder = path.abspath(path.dirname(__file__))

# Get the long description from the README and HISTORY files
with open(path.join(this_folder, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()
with open(path.join(this_folder, 'HISTORY.md'), encoding='utf-8') as f:
    history = f.read()

# Get code version from __init__.py (see https://github.com/dfm/emcee/blob/master/setup.py)
vre = re.compile("__version__ = \"(.*?)\"")
m = open(path.join(this_folder, "pygaia", "__init__.py")).read()
code_version = vre.findall(m)[0]

setup(
    name='PyGaia',

    version=code_version,

    description='Basic Gaia data simulation, manipulation, and analysis toolkit',
    long_description=long_description + "\n\n"
                     + '# Changelog\n'
                     + '\n\n'
                     + history,
    long_description_content_type="text/markdown",

    url='https://github.com/agabrown/PyGaia',

    author='Anthony Brown',
    author_email='brown@strw.leidenuniv.nl',

    license='LGPLv3+',

    # packages=['pygaia', 'pygaia.errors', 'pygaia.astrometry', 'pygaia.photometry', 'pygaia.plot',
    #  'pygaia.tests'],
    packages=find_packages(exclude=['contrib', 'docs', 'tests']),
    package_data={'': ['LICENSE', 'AUTHORS.md', 'HISTORY.md', 'INSTALL.md', 'MANIFEST.in'],
                  'pygaia': ['data/*.txt']},
    include_package_data=True,
    install_requires=['numpy', 'scipy'],

    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU Lesser General Public License v3 or later (LGPLv3+)',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Astronomy',
    ],

    entry_points={},
)
