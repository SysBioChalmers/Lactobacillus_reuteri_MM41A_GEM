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
import pandas as pd
import cobra
import My_def
import re

os.chdir('../../ComplementaryData/Step_02_DraftModels/Template/template_models/')
iBT721 = cobra.io.load_json_model('iBT721_standlized.json')
iNF517 = cobra.io.load_json_model('iNF517_standlized.json')
iML1515 = cobra.io.load_json_model('iML1515_standlized.json')

Lreu_metacyc = cobra.io.read_sbml_model('../../RAVEN/Lreu_ra_me.xml')
Lreu_kegg = cobra.io.read_sbml_model('../../RAVEN/Lreu_ra_ke.xml')
Lreu_caeveme = cobra.io.load_json_model('../../CarveMe/Lreu_ca_gp.json')
Lreu_seed = cobra.io.read_sbml_model('../../SEED/Lreu_se.xml')

Lreu_metacyc.id = 'Lreu_metacyc'
Lreu_kegg.id = 'Lreu_kegg'
Lreu_caeveme.id = 'Lreu_caeveme'
Lreu_seed.id = 'Lreu_seed'


#%%

all_map = pd.DataFrame()
#all_map = pd.read_csv('all_map.csv', sep = '\t')

for model in [Lreu_seed]:#,iNF517,iML1515,Lreu_metacyc,Lreu_keegiBT721,iNF517,iML1515,Lreu_metacyc,Lreu_kegg,Lreu_caeveme,

    print(model.id)
    if model.id in ['iBT721','iNF517','iML1515','Lreu_caeveme']:
        fromdb = 'bigg'

    elif  model.id =='Lreu_metacyc':
        fromdb = 'metacyc'
    elif model.id =='Lreu_kegg':
        fromdb = 'kegg'
    elif model.id =='Lreu_seed':
        fromdb = 'seed'

    realist = [i.id for i in model.reactions]

    metlist = [i.id for i in model.metabolites]

    trimmed_metlist = [re.sub('_..?$', '', i) for i in metlist]

    if  model.id =='Lreu_metacyc':
        trimmed_realist = [re.sub(r'__45__', r'-', i) for i in realist]
        trimmed_metlist = [re.sub(r'__45__', r'-', i) for i in metlist]
        trimmed_realist = [re.sub(r'__46__', r'.', i) for i in trimmed_realist]
        trimmed_metlist = [re.sub(r'__46__', r'.', i) for i in trimmed_metlist]
        trimmed_realist = [re.sub(r'__43__', r'+', i) for i in trimmed_realist]
        trimmed_metlist = [re.sub(r'__43__', r'+', i) for i in trimmed_metlist]
    elif model.id =='Lreu_seed':
        trimmed_realist = [re.sub(r'_[cpeb](0?)$', '', i) for i in realist]
    else:
        trimmed_realist = realist

    # %% rea
    metacyc_targetlist, MNX_IDlist = My_def.mapIDsViaMNXref('rxns',trimmed_realist,fromdb,'metacyc')
    bigg_targetlist, MNX_IDlist = My_def.mapIDsViaMNXref('rxns',trimmed_realist,fromdb,'bigg')
    kegg_targetlist, MNX_IDlist = My_def.mapIDsViaMNXref('rxns',trimmed_realist,fromdb,'kegg')
    seed_targetlist, MNX_IDlist = My_def.mapIDsViaMNXref('rxns',trimmed_realist,fromdb,'seed')

    metacyc_targetlist = [ str(i) for i  in metacyc_targetlist ]
    kegg_targetlist = [ str(i) for i  in kegg_targetlist ]
    seed_targetlist = [ str(i) for i  in seed_targetlist ]
    MNX_IDlist = [ str(i) for i  in MNX_IDlist ]


    rea_map= pd.DataFrame(list(zip( realist, bigg_targetlist ,MNX_IDlist ,metacyc_targetlist ,kegg_targetlist ,seed_targetlist)),
                     columns=['id_in_tp','bigg','metnetx','metacyc','kegg','seed'])
    rea_map['type'] = 'rea'

    # %% mets

    metacyc_targetlist, MNX_IDlist = My_def.mapIDsViaMNXref('mets',trimmed_metlist,fromdb,'metacyc')
    bigg_targetlist, MNX_IDlist = My_def.mapIDsViaMNXref('mets',trimmed_metlist,fromdb,'bigg')
    kegg_targetlist, MNX_IDlist = My_def.mapIDsViaMNXref('mets',trimmed_metlist,fromdb,'kegg')
    seed_targetlist, MNX_IDlist = My_def.mapIDsViaMNXref('mets',trimmed_metlist,fromdb,'seed')


    met_map= pd.DataFrame(list(zip( metlist, bigg_targetlist ,MNX_IDlist ,metacyc_targetlist ,kegg_targetlist ,seed_targetlist)),
                     columns=['id_in_tp','bigg','metnetx','metacyc','kegg','seed'])
    met_map['type'] = 'met'


    model_map = met_map.append(rea_map,ignore_index = True)

    model_map['model'] = model.id

    model_map.to_csv(model.id+'_map.csv',sep = '\t',index=False)

    all_map = all_map.append(model_map)

    # all_map.loc[all_map.model.isin(model_map.model)&
    #             all_map.type.isin(model_map.type)&
    #             all_map.id_in_tp.isin(model_map.id_in_tp),all_map.columns] = model_map[model_map.columns]

    #all_map = all_map.update(model_map, how = 'inner',on= ['id_in_tp','type','model'],)

    locals()[model.id+'_map'] = model_map


all_map.to_csv('all_map.csv', sep = '\t',index=False)



