#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 2019-06-14

"""add_metacyc_rea.py
:description : script
:param : 
:returns: 
:rtype: 
"""
# note: not finissed

import cobra
import os
import pandas as pd

os.chdir('../../ComplementaryData/Step_03_Compare_Refine/')

#%% load data
iNF517 = cobra.io.load_json_model('iNF517.json')
iML1515 = cobra.io.load_json_model('iML1515.json')
Lreu_merged = cobra.io.load_json_model('../Step_03_Compare_Refine/Lreu_merged_gapfiled.json')

Lreu_merged_initial = cobra.io.load_json_model('Lreu_merged.json')
Lre_metacyc = cobra.io.read_sbml_model('/Users/lhao/Box Sync/Projects/Project_Lreuteri/Lactobacillus_reuteri_MM41A_GEM/ComplementaryData/Step_02_DraftModels/RAVEN/Lreu_ra_me.xml')

universial_report2 = pd.read_csv('universial_report2.csv',sep = '\t')
metacyc_report = universial_report2[universial_report2['model'] == 'Lreu_me']




