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

from distutils.core import setup, Extension
import numpy

# the c extension module
extension_mod = Extension("pymef3", ["pymef3.c"])

setup(name = "pymef3", ext_modules=[extension_mod],
      include_dirs=[numpy.get_include()]) # This line needed for MSEL (+ the import at the beginning)