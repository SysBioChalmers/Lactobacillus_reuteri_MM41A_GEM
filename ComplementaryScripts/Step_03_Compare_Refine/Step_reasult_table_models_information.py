#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 2019-06-06

"""branchwork_model_mapping.py
:description : script
:param : 
:returns: 
:rtype: 
"""

import os
import numpy as np
import pandas as pd
import cobra

os.chdir('../../ComplementaryData')
# %% os.chdir('ComplementaryData/Step_03_Compare_Refine/')
print('----- loading data -----')
iHL622 = cobra.io.load_json_model('../ModelFiles/iHL622.json')
iNF517 = cobra.io.read_sbml_model('Initial_data/template_models/iNF517.xml')
LbReueri = cobra.io.read_sbml_model('Initial_data/template_models/Lreuteri_530.xml')
iBT721 = cobra.io.read_sbml_model('Initial_data/template_models/iBT721.xml')
iML1515 = cobra.io.read_sbml_model('Initial_data/template_models/iML1515.xml')

model_list = [iHL622, iNF517, LbReueri, iBT721, iML1515]

# %%
info_dic = {'Model': ['iHL622', 'iNF517', 'LbReuteri', 'iBT721', 'iML1515'],
            'Genes': [2019,2339,1943,3063,4243],
            'Included': [],
            'Included_persent': [],
            'Reactions': [],
            'With_GPR': [],
            'With_GPR_persent': [],
            'Common_with_iHL622': [],
            'Internal': [],
            'Transport': [],
            'Exchange': [],
            'Metabolites': [],
            'Unique': [],
            'Cytosol': [],
            'Others': []}
iHL622_rea_set = set([i.id for i in iHL622.reactions])

for model in model_list:
    info_dic['Included'].append(len(model.genes))
    info_dic['Reactions'].append(len(model.reactions))
    With_GPR = 0
    Internal = 0
    Transport = 0
    Exchange = 0
    if model.id == 'MODEL1507180045':
        for met_i in model.metabolites:
            met_i_temp = met_i.id.replace('LSQBKT','')
            met_i_temp = met_i_temp.replace('_RSQBKT','')
            met_i.id = met_i_temp
    rea_set = set([])
    for rea_i in model.reactions:
        rea_set.add(rea_i.id)
        if rea_i.gene_reaction_rule != '':
            With_GPR += 1
        if model.id == 'MODEL1507180045':
            if '_e' in rea_i.reaction and '_c' in rea_i.reaction:
                Transport += 1
        else:
            compartment_set = set([])
            for met_i in rea_i.metabolites.keys():
                compartment_set.add(met_i.compartment)
            if len(compartment_set) > 1:
                Transport += 1

    Exchange = len(model.exchanges)
    Internal = len(model.reactions) - Transport - Exchange
    info_dic['With_GPR'].append(With_GPR)
    info_dic['Common_with_iHL622'].append(len(iHL622_rea_set & rea_set))
    info_dic['Internal'].append(Internal)
    info_dic['Transport'].append(Transport)
    info_dic['Exchange'].append(Exchange)

    Metabolites = len(model.metabolites)
    Cytosol = len([i for i in model.metabolites if i.id.endswith('_c')])
    Others = Metabolites-Cytosol
    info_dic['Metabolites'].append(Metabolites)
    info_dic['Cytosol'].append(Cytosol)
    info_dic['Others'].append(Others)
    # info_dic['Unique'].append(len(set([i.id.split('_')[0] for i in model.metabolites ])))
    info_dic['Unique'] = info_dic['Cytosol']
info_dic['Included_persent'] = np.array(info_dic['Included'])/np.array(info_dic['Genes'])*100
info_dic['With_GPR_persent'] = np.array(info_dic['With_GPR'])/np.array(info_dic['Reactions'])*100


info_table = pd.DataFrame(info_dic).T
info_table.columns = info_table.iloc[0]
info_table = info_table.drop(info_table.index[0])
print(info_table)