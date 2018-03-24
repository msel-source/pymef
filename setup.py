#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  2 12:49:21 2017

setup.py file for pymef3 library

Ing.,Mgr. (MSc.) Jan Cimbálník
Biomedical engineering
International Clinical Research Center
St. Anne's University Hospital in Brno
Czech Republic
&
Mayo systems electrophysiology lab
Mayo Clinic
200 1st St SW
Rochester, MN
United States
"""

from setuptools import setup, Extension
import numpy
import sysconfig

# the c extension module
mef_file_ext = Extension("pymef.mef_file.pymef3_file",
						 ["pymef/mef_file/pymef3_file.c"],
						 extra_compile_args = ['-O3'])

setup(name = "pymef",
	version='0.1.2',
	description='Wrapper for MEF (multiscale electrophysiology format)',
      #url='http://github.com/storborg/funniest',
      author='Jan Cimbalnik',
      author_email='jan.cimbalnik@fnusa.cz, jan.cimbalnik@mayo.edu',
      #license='MIT',
      platforms = ['linux'],
      keywords='MEF Mayo electrophysiology',
      install_requires=['numpy'],
      zip_safe=False,
      classifiers=[#'License :: OSI Approved :: MIT License',
                   'Operating System :: MacOS',
                   'Operating System :: Microsoft :: Windows',
                   'Operating System :: POSIX :: Linux',
                   'Programming Language :: Python :: 2',
                   'Programming Language :: Python :: 3',
                   'Development Status :: 1 - Production/Unstable',
                   'Topic :: Scientific/Engineering'],
      packages = ["pymef","pymef.mef_file"],
      ext_modules=[mef_file_ext],
      include_dirs=[numpy.get_include()]) # This line needed for MSEL (+ the import at the beginning)