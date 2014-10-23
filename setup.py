#coding: utf8

"""
Setup script for minitti.
"""

from glob import glob


try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup


setup(name='minitti',
      version='0.0.1',
      author="TJ Lane",
      author_email="tjlane@stanford.edu",
      description='Gas pew pew',
      packages=["minitti"],
      package_dir={"minitti": "minitti"},
      scripts=[s for s in glob('scripts/*') if not s.endswith('__.py')],
      test_suite="test")
