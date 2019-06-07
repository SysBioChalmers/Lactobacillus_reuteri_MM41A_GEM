#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 2019-06-06

"""temp2universial_model.py
:description : script
:param : 
:returns: 
:rtype: 
"""

import os
import pandas as pd
import cobra
import My_def


os.chdir('../../ComplementaryData/Step_02_DraftModels/Template/template_models/')

all_map = pd.read_csv('all_map.csv', sep = '\t')
iBT721 = cobra.io.load_json_model('iBT721_standlized.json')
iNF517 = cobra.io.load_json_model('iNF517_standlized.json')
iML1515 = cobra.io.load_json_model('iML1515_standlized.json')



Lreu_metacyc = cobra.io.read_sbml_model('../../RAVEN/Lreu_ra_me.xml')
Lreu_keeg = cobra.io.read_sbml_model('../../RAVEN/Lreu_ra_ke.xml')

iBT721_initial_report = pd.read_csv('iBT721_initial_report.csv', sep='\t')
iNF517_initial_report = pd.read_csv('iNF517_initial_report.csv', sep='\t')
iML1515_initial_report = pd.read_csv('iML1515_initial_report.csv', sep='\t')

iBT721_initial_report['model'] = 'iBT721'
iML1515_initial_report['model'] = 'iML1515'
iNF517_initial_report['model'] = 'iNF517'

case = 'other'
# %%

if case =='first':
    universial_report = iBT721_initial_report.append([iML1515_initial_report,iNF517_initial_report])
    universial_report = universial_report.reset_index()

    universial_report = universial_report.fillna('')



    #mets:
    # id name Formula charge annotation
    universial_report['name'] = ''
    universial_report['formula'] = ''
    universial_report['charge'] = ''
    universial_report['annotation'] = ''

    #reas:
    # id name rea euq bounds mets gpr notes

    universial_report['metabolites'] = ''
    universial_report['reactions'] = ''
    universial_report['bounds'] = ''
    universial_report['gpr'] = ''


    columns = ['model', 'type', 'id_in_tp', 'new_id', 'id_in_bigg',
               'name','formula','charge',
               'metabolites', 'reactions','bounds', 'gpr',
               'descripation', 'old_bigg_ids', 'feature_tp', 'feature_bigg','diff', 'notes',  'annotation']


    universial_report = universial_report[columns]
    universial_report.to_csv('universalreport.csv',sep = '\t',index=False)



    # %%
    for index ,row in universial_report.iterrows():
        row = row.copy()
        if row['model'] =='iBT721':
            model = iBT721
        elif row['model'] =='iNF517':
            model = iNF517
        elif row['model'] =='iML1515':
            model = iML1515

        if row['new_id'] == '':
            id = row['id_in_tp']
        else:
            id = row['new_id']

        if row['type'] =='met':
            met = model.metabolites.get_by_id(id)

            row['name'] = met.name
            row['formula'] =  met.formula
            row['charge'] = met.charge
            row['annotation'] = str(met.annotation)

        elif row['type'] =='rea':
            rea = model.reactions.get_by_id(id)

            row['name'] = rea.name
            row['metabolites'] = My_def.model_report.str_rea_metabolites(rea.metabolites)
            row['reactions'] = rea.reaction
            row['bounds'] = str(rea.bounds)
            row['gpr'] = rea.gene_reaction_rule
            row['annotation'] = str(rea.annotation)

        universial_report.iloc[index] = row

    universial_report = universial_report.sort_values(by = ['model','type','id_in_tp'])
    universial_report.to_csv('universal_report.csv',sep = '\t',index=False)

else:
    universial_report = pd.read_csv('universal_report.csv', sep = '\t')

    universial_report2 = pd.merge(universial_report,all_map,how = 'outer',on = ['model','type','id_in_tp'])
    universial_report2 = universial_report2.reset_index()
    universial_report2 = universial_report2.fillna('')
    #universial_report2.to_csv('universial_report2.csv',sep = '\t',index=False)




