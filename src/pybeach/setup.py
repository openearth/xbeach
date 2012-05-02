from setuptools import setup, find_packages
import sys, os

version = '0.1'

setup(name='XBeach',
      version=version,
      description="XBeach is a two-dimensional model for wave propagation, long waves and mean flow, sediment transport and morphological changes of the nearshore area, beaches, dunes and backbarrier during storms.",
      long_description="""\
XBeach is a two-dimensional model for wave propagation, long waves and mean flow, sediment transport and morphological changes of the nearshore area, beaches, dunes and backbarrier during storms. It is a public-domain model that has been developed with funding and support by the US Army Corps of Engineers, by a consortium of UNESCO-IHE, Deltares (Delft Hydraulics), Delft University of Technology and the University of Miami.""",
      classifiers=[
          "Topic :: Scientific/Engineering",
          "Topic :: Scientific/Engineering :: Physics",
          "Intended Audience :: Science/Research",
          "Development Status :: 4 - Beta",
          "License :: OSI Approved :: GNU General Public License v3 (GPLv3)"
          ], # Get strings from http://pypi.python.org/pypi?%3Aaction=list_classifiers
      keywords='xbeach numerical beach coast simulation physics',
      author='Fedor Baart',
      author_email='fedor.baart@deltares.nl',
      url='http://www.xbeach.org',
      license='GPLv3',
      packages=find_packages(exclude=['ez_setup', 'examples', 'tests']),
      include_package_data=True,
      zip_safe=False,
      install_requires=[
          "numpy",
          "matplotlib",
          "nose",
          "teamcity-nose" # needed for Deltares CI.
          # -*- Extra requirements: -*-
      ],
      entry_points="""
      # -*- Entry points: -*-
      """,
      )
