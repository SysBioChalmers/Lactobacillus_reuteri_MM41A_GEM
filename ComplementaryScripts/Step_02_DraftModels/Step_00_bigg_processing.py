#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 2019-05-17

"""Step_00_bigg_processing.py
:description : script to process bigg model to easy to load
:param : bigg database
:returns: new files that easy to load
:rtype: file
"""

import os
import pickle
import cobra
import pandas as pd


def str_rea_metabolites(rea_metabolites):
    '''

    :param rea_metabolites: reaction.metabolites a dictionary
    :return: a str of rea_metabolites
    '''
    temp_dic = {}
    for k, v in rea_metabolites.items():
        temp_dic[k.id] = v
    str_mets = str(temp_dic)
    return str_mets


os.chdir('../../ComplementaryData/bigg_database/')
# %% < read bigg model>

case = 'other' # 'first' or 'other'
if case =='first':
    # pockle bigg model
    bigg_model = cobra.io.load_json_model('universal_model.json')
    with open('universal_model.pickle', 'wb') as f:
       pickle.dump(bigg_model, f)
else:
    with open('universal_model.pickle', 'rb') as f:
        bigg_model = pickle.load(f)

# %% < build met report>
bigg_met_df1 = pd.read_csv('bigg_database/bigg_models_metabolites.txt', sep='\t', usecols=[0, 1, 2, 5])
bigg_met_df1.to_csv('bigg_met_df.csv',sep = '\t',index=False)

# %% < build met report >
# the section to build all bigg_rea_df (    bigg_rea_df.to_csv('../../../bigg_database/bigg_models_reactions_all.txt'))
# becaues the reaction_strings in file are not irregular
# bigg_rea_df from file
#   option 1
bigg_rea_df1 = pd.read_csv('bigg_models_reactions.txt', sep='\t', usecols=[0, 2, 5])

#   option 2 bigg_rea_df from model
columns = ['bigg_id', 'equaction_string', 'mets', 'bounds', 'original_bigg_ids']
bigg_rea_df2 = pd.DataFrame(columns=columns)

for rea in bigg_model.reactions:
    rowlist = []

    rowlist.append(rea.id)
    rowlist.append(rea.reaction)
    rowlist.append(str_rea_metabolites(rea.metabolites))
    rowlist.append(rea.bounds)

    if 'original_bigg_ids' in rea.notes.keys():
        original_bigg_ids = rea.notes['original_bigg_ids']
    else:
        original_bigg_ids = []
    rowlist.append(original_bigg_ids)

    bigg_rea_df2.loc[len(bigg_rea_df2)] = rowlist
#bigg_rea_df2.to_csv('bigg_models_reactions_from_model.txt', index=False)
bigg_rea_df = pd.merge(bigg_rea_df1, bigg_rea_df2, how='outer', on='bigg_id')
bigg_rea_df.to_csv('bigg_rea_df.csv', sep = '\t', index=False)



