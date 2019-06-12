#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 2019-06-12

"""io_outmetatooldata.py
:description : script
:param : 
:returns: 
:rtype: 
"""
import cobra.test
model = cobra.test.create_test_model("textbook")
filename = '/Users/lhao/Documents/MATLAB/aumic-master/aumic/text.dat'
f  = open(filename,'w')




f.write('#test\n')

ENZREV = []
ENZIRREV = []
for i in model.reactions:
    if i.reversibility:
        ENZREV.append(i.id)
    else:
        ENZIRREV.append(i.id)

METINT = []
METEXT = []

for i in model.metabolites:
    if i.compartment == 'c':
        METINT.append(i.id)
    else:
        METEXT.append(i.id)

f.write('-ENZREV\n')
f.write(' '.join(ENZREV))

f.write('-ENZIRREV\n')
f.write(' '.join(ENZIRREV))

f.write('-METINT\n')
f.write(' '.join(METINT))

f.write('-METEXT\n')
f.write(' '.join(METEXT))











