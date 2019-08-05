#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    use_scm_version=True,
    setup_requires=['setuptools_scm'],
    install_requires=['setuptools_scm'],
    tests_require=['pytest', 'flake8'],
    name='runControl',
    author='Maryam Zaheri',
    author_email='lastname@gmail.com',
    description='Run control analysis for NGS run with HIV sample.',
    packages=find_packages('src'),
    package_dir={'': 'src'},
    entry_points={
        'console_scripts': ['runControl = runControl.cli:main']
    },
    license='MIT',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/medvir/runControl",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
