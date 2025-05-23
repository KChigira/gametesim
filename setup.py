#!/usr/bin/env python
from setuptools import setup, find_packages
from gametesim.__init__ import __version__

setup(name='gametesim',
      version=__version__,
      description='Simulating recombination by crossing in plant breeding.',
      author='Koki Chigira',
      author_email='chigirak@g.ecc.u-tokyo.ac.jp',
      url='https://github.com/KChigira/gametesim/',
      license='MIT',
      packages=find_packages(),
      install_requires=[
      ],
      entry_points={'console_scripts': [
            'parents = gametesim.parents:main',
            'cross = gametesim.cross:main',
            'selfing = gametesim.selfing:main',
            'backcross = gametesim.backcross:main',
            'randcross = gametesim.randcross:main',
            'genostat = gametesim.genostat:main',
            'genovisual = gametesim.genovisual:main',
            'genomap = gametesim.genomap:main',
            ]
      }
    )
