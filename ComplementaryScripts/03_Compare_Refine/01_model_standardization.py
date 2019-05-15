#!/usr/bin/env python
#-*- coding: utf-8 -*-
# Created by lhao at 2019-04-18

import os
import cobra
import re
import My_def
import pandas as pd

os.chdir('/Users/lhao/Box Sync/Project_Lreuteri/03_Compare_Refine/')


namelist = ['Lreu_se_rast.xml','Lreu_se_pro.xml','Lreu_se_dna.xml']
for i in namelist:
    model = cobra.io.read_sbml_model('initial_models/' + i)
    for reas in model.reactions:
        if str(reas.lower_bound)  == '-inf':
            reas.lower_bound = -1000
            reas.upper_bound = 1000
    cobra.io.save_json_model(model,'initial_models/'+ i.split('.')[0]+'.json')



if False:

    Lreu_ra_te = cobra.io.read_sbml_model('initial_models/Lreu_ra_te.xml')
    Lreu_py_te = cobra.io.read_sbml_model('initial_models/Lreu_py_te.xml')
    iNF517 = cobra.io.load_json_model('iNF517.json')

    def add_inf(submodel,iNF517):

        for i in submodel.reactions:
            initial_i = iNF517.reactions.get_by_id(i.id)
            i.notes.update(initial_i.notes)
            i.annotation.update(initial_i.annotation)

        for i in submodel.metabolites:
            initial_i = iNF517.metabolites.get_by_id(i.id)
            i.notes.update(initial_i.notes)
            i.annotation.update(initial_i.annotation)
        return

    add_inf(Lreu_ra_te,iNF517)
    add_inf(Lreu_py_te,iNF517)

    cobra.io.save_json_model(Lreu_ra_te,'Lreu_ra_te.json')
    cobra.io.save_json_model(Lreu_py_te,'Lreu_py_te.json')

    print(Lreu_ra_te.reactions[10].notes)
    print(Lreu_py_te.reactions[10].notes)

    Lreu_ra_te = cobra.io.load_json_model('Lreu_ra_te.json')
    Lreu_py_te = cobra.io.load_json_model('Lreu_py_te.json')




