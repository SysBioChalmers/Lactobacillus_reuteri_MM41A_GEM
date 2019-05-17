#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by lhao at 2019-05-17

'''
input: L.reuteri protein sequence
output: draft model

'''

import os
import cobra
import My_def

os.chdir('../../ComplementaryData/02_DraftModels/')

#Gram positive
os.system('carve Lreuteri_biogaia_v03.faa --cobra -u grampos -o CarveMe/Lreu_ca_gp.xml');

model = cobra.io.read_sbml_model('CarveMe/Lreu_ca_gp.xml')
model.description = 'GEM of L reuteri by CarveMe'
cobra.io.save_json_model(model,'CarveMe/Lreu_ca_gp.json')
My_def.io_outtxt(model,'CarveMe/Lreu_ca_gp.txt',True)


#all
os.system('carve Lreuteri_biogaia_v03.faa --cobra  -o CarveMe/Lreu_ca.xml');

model = cobra.io.read_sbml_model('CarveMe/Lreu_ca.xml')
model.description = 'GEM of L reuteri by CarveMe'
cobra.io.save_json_model(model,'CarveMe/Lreu_ca.json')
My_def.io_outtxt(model,'CarveMe/Lreu_ca.txt',True)
