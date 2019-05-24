#!/usr/bin/env python
#-*- coding: utf-8 -*-
# Created by lhao at 2019-04-15


import os
import cobra
import re
import My_def
import pandas as pd
import pickle


os.chdir('/Users/lhao/Box Sync/Project_Lreuteri/Step_03_Compare_Refine/')


def full_list(mnxlist,fulllist):
    gaplist = []
    for i, v in enumerate(mnxlist):
        if v =='':
            gaplist.append(fulllist[i])
            mnxlist[i] = fulllist[i]
    return mnxlist,gaplist

modelsnamelist = ['Lreu_ca.xml','Lreu_ca_gp.xml',
             'Lreu_se_rast_cobra.json','Lreu_se_pro_cobra.json','Lreu_se_dna_cobra.json',
             'Lreu_ra_ke.xml','Lreu_ra_me.xml','Lreu_ra_te.json',
             'Lreu_py_te.json',
             ]

modelsidspace =  ['bigg','bigg',
                  'seed','seed','seed',
                  'kegg','metacyc','bigg',
                  'bigg'
                ]

modelsdescribe = ['CarveMe','CarveMe ues Gram positive template model',
                  'SEED Rast','SEED pro seq','SEED DNA seq',
                  'RAVEN KEEG','RAVEN MetaCyc','RAVEN template',
                  'Python template'
                ]

modelpdlist = []
for model_index in range(len(modelsnamelist)):

    if modelsnamelist[model_index].endswith(r'.xml'):
        model = cobra.io.read_sbml_model('initial_models/'+modelsnamelist[model_index])
    elif modelsnamelist[model_index].endswith(r'.json'):
        model = cobra.io.load_json_model('initial_models/'+modelsnamelist[model_index])


    gene_set = set([i.id for i in model.genes])
    realist = [i.id for i in model.reactions]
    metlist = [i.id for i in model.metabolites]


    model.notes['idspace'] = modelsidspace[model_index]
    model.notes['describe'] = modelsdescribe[model_index]


    if model_index in [0,1]:
        pattern = r'_[cpe]$'
        trimmed_realist = realist
        trimmed_metlist = [re.sub(pattern, '', i) for i in metlist]

    elif model_index in [2,3,4]:
        pattern = r'_[cpeb](0?)$'
        trimmed_realist = [re.sub(pattern, '', i) for i in realist]
        trimmed_metlist = [re.sub(pattern, '', i) for i in metlist]

    elif  model_index ==5:
        trimmed_realist = realist
        trimmed_metlist = metlist

    elif model_index ==6:
        trimmed_realist = [re.sub(r'__45__', r'-', i) for i in realist]
        trimmed_metlist = [re.sub(r'__45__', r'-', i) for i in metlist]
        trimmed_realist = [re.sub(r'__46__', r'.', i) for i in trimmed_realist]
        trimmed_metlist = [re.sub(r'__46__', r'.', i) for i in trimmed_metlist]
        trimmed_realist = [re.sub(r'__43__', r'+', i) for i in trimmed_realist]
        trimmed_metlist = [re.sub(r'__43__', r'+', i) for i in trimmed_metlist]

    elif model_index in [7,8]:
        pattern = r'_[cpe]$'

        trimmed_realist = realist
        trimmed_metlist = [re.sub(pattern, '', i) for i in metlist]

        _, realist_mnx = My_def.mapIDsViaMNXref('rxns', trimmed_realist, model.notes['idspace'], 'metanetx')
        _, metlist_mnx = My_def.mapIDsViaMNXref('mets', trimmed_metlist, model.notes['idspace'], 'metanetx')
        print(realist_mnx.count(''), metlist_mnx.count(''))

        indexrealist = [i for i, x in enumerate(realist_mnx) if x == '']
        indexmetlist = [i for i, x in enumerate(metlist_mnx) if x == '']

        for i in indexrealist:
            trimmed_realist[i] = model.reactions[i].notes['original_bigg_ids'][0]
        for i in indexmetlist:
            temp = model.metabolites[i].notes['original_bigg_ids'][0]
            trimmed_metlist[i] = re.sub(pattern, '', temp)

        _, realist_mnx = My_def.mapIDsViaMNXref('rxns', trimmed_realist, model.notes['idspace'], 'metanetx')
        print(realist_mnx.count(''), metlist_mnx.count(''))

        indexrealist = [i for i, x in enumerate(realist_mnx) if x == '']

        for i in indexrealist:
            trimmed_realist[i] = model.reactions[i].notes['original_bigg_ids'][-1]

    realist_mnx = My_def.mapIDsViaMNXref('rxns', trimmed_realist, model.notes['idspace'], 'metanetx')[1]
    metlist_mnx = My_def.mapIDsViaMNXref('mets', trimmed_metlist, model.notes['idspace'], 'metanetx')[1]

    rea_set, reagaplist = full_list(realist_mnx, trimmed_realist)
    met_set, metgaplist = full_list(metlist_mnx, trimmed_metlist)

    for i in range(len(realist)):
        model.reactions[i].notes['metanetx_id']   = realist_mnx[i]
    for i in range(len(metlist)):
        model.metabolites[i].notes['metanetx_id'] = metlist_mnx[i]

    print(modelsnamelist[model_index].split('.')[0])
    print(len(reagaplist), len(metgaplist))
    print(reagaplist)
    print(metgaplist)

    cobra.io.save_json_model(model,modelsnamelist[model_index].split('.')[0]+'.json')

    #modelnotepd = pd.DataFrame({'trimmed_realist':trimmed_realist})

    pddata = [

        pd.DataFrame({'realist': realist}),
        pd.DataFrame({'trimmed_realist': trimmed_realist}),
        pd.DataFrame({'rea_set': list(rea_set)}),
        pd.DataFrame({'reagaplist': reagaplist}),

        pd.DataFrame({'metlist': metlist}),
        pd.DataFrame({'trimmed_metlist': trimmed_metlist}),
        pd.DataFrame({'met_set': list(met_set)}),
        pd.DataFrame({'metgaplist': metgaplist}),

        pd.DataFrame({'gene_set': list(gene_set)}),


    ]

    modelnotepd = pd.concat(pddata,axis = 1)
    modelnotepd = modelnotepd.fillna('')

    model_var_name_str = modelsnamelist[model_index].split('.')[0] + 'notepd'
    locals()[model_var_name_str] = modelnotepd

    modelpdlist.append(locals().get(model_var_name_str))


with open('data.pickle', 'wb') as f:
    pickle.dump(modelpdlist, f)

print('Done')








