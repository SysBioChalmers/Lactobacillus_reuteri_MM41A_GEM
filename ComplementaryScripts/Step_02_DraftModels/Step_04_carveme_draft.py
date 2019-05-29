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
import pandas as pd



os.chdir('../../ComplementaryData/Step_02_DraftModels/')
case = 'other' #'first' or 'other'

# %% <build>
if case =='frist':
    #Gram positive
    os.system('carve Lreuteri_biogaia_v03.faa --cobra -u grampos -o CarveMe/Lreu_ca_gp.xml');

    #all
    os.system('carve Lreuteri_biogaia_v03.faa --cobra  -o CarveMe/Lreu_ca.xml');

# %% <standstandardlization>

def CarveMe_processing(covermemodel):
    #change gene id 'G_id'
    for gen in  covermemodel.genes:
        gen.id = gen.id.replace('G_','')
    # combine met according report
    My_def.model_report.combine_met('cyst__L_c','cysth__L_c',covermemodel)

    return covermemodel

Lreu_ca_gp = cobra.io.read_sbml_model('CarveMe/Lreu_ca_gp.xml')
Lreu_ca_gp.description = 'GEM of L reuteri by CarveMe'
Lreu_ca_gp.id = 'Lreu_ca_gp'

Lreu_ca = cobra.io.read_sbml_model('CarveMe/Lreu_ca.xml')
Lreu_ca.description = 'GEM of L reuteri by CarveMe'
Lreu_ca.id = 'Lreu_ca'



bigg_rea_df = pd.read_csv('../bigg_database/bigg_rea_df.csv', sep='\t')
bigg_met_df = pd.read_csv('../bigg_database/bigg_met_df.csv', sep='\t')

Lreu_ca_standardlized, ca_report = My_def.model_report.model_report_compare_bigg(Lreu_ca, bigg_rea_df, bigg_met_df, compartment='_')
Lreu_ca_gp_standardlized, ca_gp_report = My_def.model_report.model_report_compare_bigg(Lreu_ca_gp, bigg_rea_df, bigg_met_df,
                                                                   compartment='_')
# %% <Manual change according the report>

Lreu_ca_gp_standardlized = CarveMe_processing(Lreu_ca_gp_standardlized)
Lreu_ca_standardlized = CarveMe_processing(Lreu_ca_standardlized)

cobra.io.save_json_model(Lreu_ca_standardlized, 'CarveMe/Lreu_ca.json')
cobra.io.save_json_model(Lreu_ca_gp_standardlized, 'CarveMe/Lreu_ca_gp.json')


#My_def.io_outtxt(Lreu_ca,'CarveMe/Lreu_ca.txt',True)
#My_def.io_outtxt(Lreu_ca_gp,'CarveMe/Lreu_ca_gp.txt',True)
