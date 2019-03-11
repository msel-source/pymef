#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 22:44:21 2019

Script to parse meflib.h file to extract MEF constants

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


f = open('/home/jan_cimbalnik/Dropbox/Source/C/pymef/meflib/meflib/meflib.h', 'r')

lines = f.readlines()

f.close()

constants_start = [i for i, x in enumerate(lines) if 'MEF Constants' in x][0]
constants_stop = [i for i, x in enumerate(lines) if 'MEF Macros' in x][0]

# Get Mef constant definitions
definitions = []
definitions_dict = {}
for line in lines[constants_start: constants_stop]:
    if '#define' in line:
        definition = line[:]
        # Remove define
        definition = definition.replace('#define', '')
        # Remove comments
        if '//' in definition:
            definition = definition[:definition.index('//')]
        # Skip if the whole line was commented
        if len(definition) == 0:
            continue
        # Remove tabs
        definition = definition.replace('\t', ' ')
        # Remove whitespace
        definition = definition.strip()
        
        # Skip mef globals
        if definition.startswith('MEF_GLOBALS'):
            continue
        
        # Skip file processing struct constants
        if definition.startswith('FPS'):
            continue
        
        key = definition[:definition.index(' ')]
        value = definition[definition.index(' '):].strip()
        
        try:
            int(value, 0)
            line_str = " = ".join([key, value])
            definitions_dict[key] = int(value, 0)
        except ValueError:
            if value in definitions_dict.keys():
                line_str = " = ".join([key, str(definitions_dict[value])])
            else:
                line_str = " = ".join([key, '"'+value+'"'])
    
        definitions.append(line_str+'\n')
             
                
# No write the contents into mef_constants.py
f = open('/home/jan_cimbalnik/Dropbox/Source/C/pymef/pymef/mef_constants.py',
         'w')

f.writelines(definitions)
f.close()
    

