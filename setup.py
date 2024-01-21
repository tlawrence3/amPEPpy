import re
import os
from setuptools import setup, Extension
from codecs import open
from os import path

version_file = open("amPEPpy/_version.py", "r").read()
version_match = re.match(r"^__version__ = ['\"]([^'\"]*)['\"]", version_file)
if (version_match):
    version = version_match.group(1)
else:
    raise RuntimeError("Unable to find version string in _version.py")

here = path.abspath(path.dirname(__file__))
# Get the long description from the README file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
        long_description = f.read()

setup(name = "amPEPpy",
      setup_requires=['pytest-runner'],
      tests_require=['pytest'],
      install_requires=['biopython', 'numpy', 'scikit-learn', 'pandas'],
      packages = ["amPEPpy"],
      python_requires='~=3.5',
      entry_points = {
          "console_scripts": ['ampep = amPEPpy.amPEP:main']},
      version = version,
      author="Dana L. Carper, Travis J. Lawrence, Alyssa A. Carrell, David J. Weston",
      author_email="",
      description = "amPEP ...",
      long_description=long_description,
      license='GPLv3',
      url = "",
      classifiers=['Development Status :: 4 - Beta',
                   'Environment :: Console',
                   'Intended Audience :: Science/Research',
                   'License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)',
                   'Natural Language :: English',
                   'Operating System :: MacOS :: MacOS X',
                   'Operating System :: POSIX :: Linux',
                   'Programming Language :: Python :: 3.6',
                   'Programming Language :: Python :: 3.7',
                   'Programming Language :: Python :: 3.8',
                   'Programming Language :: Python :: 3.9',
                   'Programming Language :: Python :: 3.10',
                   'Programming Language :: Python :: Implementation :: CPython',
                   'Topic :: Scientific/Engineering :: Bio-Informatics'],)
