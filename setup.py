#!/usr/bin/env python3

from setuptools import setup

with open('VERSION','r') as _f:
    __version__ = _f.read().strip()

setup(
    name='lnls',
    version=__version__,
    author='lnls-fac',
    description='LNLS utilities',
    url='https://github.com/lnls-fac/lnls',
    download_url='https://github.com/lnls-fac/lnls',
    license='MIT License',
    classifiers=[
        'Intended Audience :: Science/Research',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering'
    ],
    packages=['lnls'],
    package_data={'lnls': ['VERSION']},

    install_requires=[
        'numpy>=1.8.2',
        'matplotlib>=1.4.2'
    ]
)
