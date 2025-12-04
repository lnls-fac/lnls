#!/usr/bin/env python-sirius
"""Setup module."""

import pathlib

from setuptools import setup


def get_abs_path(relative):
    """."""
    return str(pathlib.Path(__file__).parent / relative)


with open(get_abs_path("README.md"), "r") as _f:
    _long_description = _f.read().strip()

with open(get_abs_path("VERSION"), "r") as _f:
    __version__ = _f.read().strip()

with open(get_abs_path("requirements.txt"), "r") as _f:
    _requirements = _f.read().strip().split("\n")


setup(
    name='lnls',
    version=__version__,
    author='lnls-fac',
    description='LNLS utilities',
    long_description=_long_description,
    url='https://github.com/lnls-fac/lnls',
    download_url='https://github.com/lnls-fac/lnls',
    license='MIT License',
    classifiers=[
        'Intended Audience :: Science/Research',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering'
    ],
    packages=['lnls'],
    install_requires=_requirements,
    package_data={'lnls': ['VERSION']},
    zip_safe=False,
)
