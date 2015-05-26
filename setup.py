#!/usr/bin/env python3

from setuptools import setup, find_packages

with open('VERSION','r') as _f:
    __version__ = _f.read().strip()

setup(
    name='lnls',
    version=__version__,
    description='LNLS utilities',
    url='https://github.com/lnls-fac/lnls',
    package_dir={'lnls': 'src'},
    packages=['lnls'],

    install_requires=[
        'numpy>=1.8.2',
        'matplotlib>=1.4.2'
    ]
)
