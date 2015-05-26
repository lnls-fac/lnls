#!/usr/bin/env python3

from distutils.core import setup

with open('VERSION', 'r') as _f:
    __version__ = _f.read().strip()

setup(
    name='LNLS',
    version=__version__,
    description='LNLS utilities',
    url='https://github.com/lnls-fac/lnls',
    packages=['lnls'],
)
