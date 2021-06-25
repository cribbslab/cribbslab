import sys
import os
import re
import setuptools
from setuptools import setup, find_packages, Extension

from distutils.version import LooseVersion
if LooseVersion(setuptools.__version__) < LooseVersion('1.1'):
    print("Version detected:", LooseVersion(setuptools.__version__))
    raise ImportError(
        "the scflow requires setuptools 1.1 higher")

########################################################################
########################################################################
IS_OSX = sys.platform == 'darwin'

########################################################################
########################################################################
# collect version
print(sys.path.insert(0, "cribbslab"))
import version

version = version.__version__

###############################################################
###############################################################
# Define dependencies
#
major, minor1, minor2, s, tmp = sys.version_info

if major < 3:
    raise SystemExit("""Requires Python 3 or later.""")

cribbslab_packages = find_packages()
cribbslab_package_dirs = {'cribbslab': 'cribbslab'}

##########################################################
##########################################################
# Classifiers
classifiers = """
Development Status :: 3 - Alpha
Intended Audience :: Science/Research
Intended Audience :: Developers
License :: OSI Approved
Programming Language :: Python
Topic :: Software Development
Topic :: Scientific/Engineering
Operating System :: POSIX
Operating System :: Unix
Operating System :: MacOS
"""

setup(
    # package information
    name='cribbslab',
    version=version,
    description='cribbslab : cribbslab workflows',
    author='Adam Cribbs',
    author_email='adam.cribbs@ndorms.ox.ac.uk',
    license="MIT",
    platforms=["any"],
    keywords="computational genomics",
    long_description='cribbslab : cribbslab workflows',
    classifiers=[_f for _f in classifiers.split("\n") if _f],
    url="",
    # package contents
    packages=cribbslab_packages,
    package_dir=cribbslab_package_dirs,
    include_package_data=True,
    entry_points={
        "console_scripts": ["cribbslab = cribbslab.entry:main"]
    },
    # other options
    zip_safe=False,
    test_suite="tests",
)
