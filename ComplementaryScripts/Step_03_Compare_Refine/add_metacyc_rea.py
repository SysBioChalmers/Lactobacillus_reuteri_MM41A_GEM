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
import re
import My_def
import pickle
from importlib import  reload
import sys
# Disable
def blockPrint():
    sys.stdout = open(os.devnull, 'w')
# Restore
def enablePrint():
    sys.stdout = sys.__stdout__


#reload(My_def)

os.chdir('../../ComplementaryData/Step_03_Compare_Refine/')

#%% load data

Lreu_draft_3_refined = cobra.io.load_json_model('../Step_02_DraftModels/Lreu_draft_3_refined.json')
Lreu_metset = set([i.id for i in Lreu_draft_3_refined.metabolites ])
Lreu_reaset = set([i.id for i in Lreu_draft_3_refined.reactions ])

Lreu_metacyc = cobra.io.read_sbml_model('/Users/lhao/Box Sync/Projects/Project_Lreuteri/Lactobacillus_reuteri_MM41A_GEM/ComplementaryData/Step_02_DraftModels/RAVEN/Lreu_ra_me.xml')

Lreu_metacyc.id = 'Lreu_metacyc'

all_map = pd.read_csv('../Step_02_DraftModels/all_map.csv', sep = '\t')
all_map = all_map.fillna('')

Lreu_metacyc_map = all_map[all_map.model == 'Lreu_metacyc'].copy()


iML1515_report = all_map[(all_map.model == 'iML1515')].copy()
iML1515_reaset = set(iML1515_report[(iML1515_report.type == 'rea')].id_in_tp.values)
iML1515_metset = set(iML1515_report[(iML1515_report.type == 'met')].id_in_tp.values)

iNF517_report = all_map[(all_map.model == 'iNF517')].copy()
iNF517_reaset = set(iNF517_report[(iNF517_report.type == 'rea')].id_in_tp.values)
iNF517_metset = set(iNF517_report[(iNF517_report.type == 'met')].id_in_tp.values)


with open('../bigg_database/universal_model.pickle', 'rb') as f:
        bigg_model = pickle.load(f)
bigg_reaset = set([i.id for i in bigg_model.reactions])
bigg_met_set = set([i.id for i in bigg_model.metabolites])



# %% compare id_map and rea_map

def find_one(str,modelset ,bigg_model):

    if str =='':
        return ''
    else:
        str = eval(str)
        biggids = [i for i in str if not (i.startswith('R_') or i.startswith('M_'))]

        if len(biggids) == 1:
            return biggids[0]
        else:
            a = set(biggids) & modelset
            if len(a)==1:
                return list(a)[0]
            elif len(a)>1 :
                biggids = list(a)
            for i in biggids:
                try:
                    rea = bigg_model.reactions.get_by_id(i)
                    if '_p' not in rea.reaction and '_m' not in rea.reaction:
                        return i
                except KeyError:
                    continue

            return biggids[0]


# %% <option1> reaid map, metacyc-->bigg  find the euqation
temp_rea_df = Lreu_metacyc_map[Lreu_metacyc_map.type == 'rea'].copy()


modelset = iML1515_reaset|iNF517_reaset

temp_rea_df['uniqueid'] = temp_rea_df.bigg.apply(lambda x: find_one(x,modelset,bigg_model))

temp_rea_df = temp_rea_df[temp_rea_df.uniqueid != '']

#uniquelist1 = set(temp_rea_df['uniqueid'].values)   #-Lreu_reaset
#temp_rea_df = temp_rea_df[~temp_rea_df['uniqueid'].isin(uniquelist1) ]

# %% change model:
remove_rea_list_1 = list(temp_rea_df.id_in_tp)
add_rea_list_1 = list(temp_rea_df.uniqueid)

Lrue_metacyc_bigg_rearep = Lreu_metacyc.copy()


for row in temp_rea_df.itertuples():
    try:
        rea_a = bigg_model.reactions.get_by_id(row.uniqueid)
    except KeyError:
        print(row)
        continue
    rea_r = Lrue_metacyc_bigg_rearep.reactions.get_by_id(row.id_in_tp)
    rea_a.gene_reaction_rule = rea_r.gene_reaction_rule
    rea_a.notes['from'] = ['metacyc']
    rea_a.notes['metacyc_id'] = row.id_in_tp
    Lrue_metacyc_bigg_rearep.add_reaction(rea_a)
    rea_r.remove_from_model()



# %% <option2> metaid map, metacyc-->bigg  get the euqation
temp_met_df = Lreu_metacyc_map[Lreu_metacyc_map.type == 'met'].copy()

modelset = iML1515_metset|iNF517_metset

temp_met_df['uniqueid'] = temp_met_df.bigg.apply(lambda x: find_one(x,modelset,bigg_model))

temp_met_df = temp_met_df[temp_met_df.uniqueid != '']

Lrue_metacyc_bigg_metarep = Lreu_metacyc.copy()

remove_met_list = list(temp_met_df.id_in_tp.values)
add_met_list = list(temp_met_df.uniqueid.values)

remove_met_set = set(remove_met_list)
add_met_set = set(add_met_list)

remove_rea_list_2 = []

for rea in Lrue_metacyc_bigg_metarep.reactions:
    rea_metset = set(i.id for i in rea.metabolites.keys())
    #print(rea_metset)

    rea_intersection = rea_metset & remove_met_set
    #print(rea_intersection)

    if rea_intersection == set():
        continue
    elif rea_intersection == rea_metset:
        remove_rea_list_2.append(rea.id)
    for i in rea_intersection:
        add_met =  add_met_list[remove_met_list.index(i)]+'_c'
        if add_met in bigg_met_set:
            rea.reaction = re.sub(r'\b%s\b'%i , add_met_list[remove_met_list.index(i)]+'_c', rea.reaction)
        else:
            print(add_met)

# %% add reaction intomode

Lreu_draft_3_add_metacyc = Lreu_draft_3_refined.copy()
dup_set_checked = {'TRANS__45__RXN0__45__617',
             'NTP1',
             'TRANS__45__RXN__45__220',
             'TRANS__45__RXN__45__366',
             'TRANS__45__RXN0__45__510',
             'TRANS__45__RXN__45__290',
             '3__46__6__46__3__46__3__45__RXN',
             'TRANS__45__RXN__45__237',
             'RXN__45__19779',
             'PHOSPHOKETOLASE__45__RXN',
             'RXN__45__16819',
             'RXN0__45__5213'}

#list1:
final_add_list1 = []
for i in set(add_rea_list_1) - Lreu_reaset - dup_set_checked:
    rea = Lrue_metacyc_bigg_rearep.reactions.get_by_id(i)
    if '_m' not in rea.reaction and '_p' not in  rea.reaction:
        final_add_list1.append(rea.id)
        print(rea)
    # rea.notes['from'] = ['metacyc']
    # rea_a.notes['metacyc_id'] = row.id_in_tp
    Lreu_draft_3_add_metacyc.add_reaction(rea)

#list2:
final_add_list2 = []
for i in set(remove_rea_list_2) - set(remove_rea_list_1) - dup_set_checked:
    rea = Lrue_metacyc_bigg_metarep.reactions.get_by_id(i)
    final_add_list2.append(rea.id)
    print(rea)
    rea.notes['from'] = ['metacyc']
    Lreu_draft_3_add_metacyc.add_reaction(rea)


# %% check duplication
# get dup_set_checked and avoide to add it to model
check_df = My_def.model_refine.remove_dup_rea(Lreu_draft_3_add_metacyc, remove = False ,skip_met = [])
#print(check_df)


check_id = list(check_df['id'])
print('Duplicate reactions: ')
for i in list(check_id):
    rea = Lreu_draft_3_add_metacyc.reactions.get_by_id(i)
    print(rea,rea.bounds,rea.gene_reaction_rule,rea.notes)
# dup_set_checked = {'TRANS__45__RXN0__45__617',
#              'NTP1',
#              'TRANS__45__RXN__45__220',
#              'TRANS__45__RXN__45__366',
#              'TRANS__45__RXN0__45__510',
#              'TRANS__45__RXN__45__290',
#              '3__46__6__46__3__46__3__45__RXN',
#              'TRANS__45__RXN__45__237',
#              'RXN__45__19779',
#              'PHOSPHOKETOLASE__45__RXN',
#              'RXN__45__16819',
#              'RXN0__45__5213'}


# %% <save >
cobra.io.save_json_model(Lreu_draft_3_add_metacyc,'Lreu_draft_3_add_metacyc.json',sort='True')


print('=====  Done =====')










