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
mef_file_ext = Extension("pymef.mef_file.pymef3_file", ["pymef/mef_file/pymef3_file.c"])

setup(name = "pymef",
      packages = ["pymef","pymef.mef_file"],
      ext_modules=[mef_file_ext],
      include_dirs=[numpy.get_include()]) # This line needed for MSEL (+ the import at the beginning)