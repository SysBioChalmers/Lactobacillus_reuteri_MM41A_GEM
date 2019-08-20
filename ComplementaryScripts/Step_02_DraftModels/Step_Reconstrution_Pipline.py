#!/usr/bin/env python
#-*- coding: utf-8 -*-
# Created by lhao at 2019-02-27


"""pipeline.py
:description : script to renconstract GEM for L reuteri
:param : template models and seqs (standardizated) iNF517,iML1515,iBT721
:returns: draft model : Lreu_draft_3_refined
:rtype:
"""

import os
import pandas as pd
import re
import My_def
import cobra
from importlib import  reload
reload(My_def)

pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)

def get_Lreu_from_tp(qseq_file,sseq_file,tp_model):

    # <step blast(used dimond) >
    blast_result_s_in_q, blast_result_q_in_s = My_def.seq_ana.blastp_pairwise(qseq_file, sseq_file, out_dir='blast/')
    # <step select blast results: cut off bitscore=100, ppos=45 best(BBHs)>
    blast_result_df = My_def.seq_ana.select_blast(blast_result_s_in_q, blast_result_q_in_s, best_match=True, evalue=10 ** -10, pident=0, length=0,
                                          bitscore=100, ppos=45, qcovs=0)
    #blast_result_df.to_csv('result.csv')

    os.system('rm ' + blast_result_s_in_q )
    os.system('rm ' + blast_result_q_in_s)
    # <step get draft model from template>
    # note: only extracted the matched reactions
    Lreu_py_tp = My_def.get_draft_from_template(tp_model, blast_result_df,remove_missing_genes = False)
    return Lreu_py_tp


def add_inf517_etb_rea(Lreu_draft_2,iNF517):
    reaset = set([i.id for i in Lreu_draft_2.reactions])
    metset = set([i.id for i in Lreu_draft_2.metabolites])
    for rea in iNF517.reactions:
        if ('EX_' in rea.id):
            # case exchange reactions
            if rea.id not in reaset:
                rea.notes['from'] = ['iNF517','EX']
                Lreu_draft_2.add_reaction(rea)
                reaset.add(rea.id)
        elif ('_c' in rea.reaction) and ('_e' in rea.reaction):
            # case transport
            if rea.id not in reaset :
                rea.notes['from'] = ['iNF517','transport']
                Lreu_draft_2.add_reaction(rea)
                reaset.add(rea.id)
        elif ('_LRE' in rea.id) or ('LRE_c' in rea.reaction):
            # case biomass
            #print(rea)
            if rea.functional:
                if rea.id not in reaset:
                    rea.notes['from'] = ['iNF517','biomass']
                    Lreu_draft_2.add_reaction(rea)
                    reaset.add(rea.id)
        elif rea.id in ['ATPM'] and rea.id not in reaset:
            #case ATPM
            rea.notes['from'] = ['iNF517','atp']
            Lreu_draft_2.add_reaction(rea)
            reaset.add(rea.id)

    return Lreu_draft_2


# %% <general step: load data: seqs, template models >

print('----- loading data -----')
os.chdir('../../ComplementaryData/Step_02_DraftModels/Template/')
os.system('cp ../../Step_01_Sequences_analysis/Lreuteri_biogaia_v03/Lreuteri_biogaia_v03.faa ../Lreuteri_biogaia_v03_2.faa')

Lreu_seq = '../Lreuteri_biogaia_v03_2.faa'

iNF517_seq = 'template_seqs/iNF517.faa'
iBT721_seq = 'template_seqs/iBT721.faa'
iML1515_seq = 'template_seqs/iML1515.faa'
Lreuteri_530_seq = 'template_seqs/Lreuteri_530.faa'

iNF517 = cobra.io.load_json_model('template_models/iNF517_standlized.json')
iBT721 = cobra.io.load_json_model('template_models/iBT721_standlized.json')
iML1515 = cobra.io.load_json_model('template_models/iML1515_standlized.json')
Lreuteri_530 = cobra.io.load_json_model('template_models/Lreuteri_530_standlized.json')


# %% <general step: blast and get draft models. steps details could be find in in get_Lreu_from_tp function (above)>

print('----- blasting and geting draft models  -----')
Lreu_from_iNF517 = get_Lreu_from_tp(Lreu_seq,iNF517_seq,iNF517)
Lreu_from_iBT721 = get_Lreu_from_tp(Lreu_seq,iBT721_seq,iBT721)
Lreu_from_iML1515 = get_Lreu_from_tp(Lreu_seq,iML1515_seq,iML1515)
Lreu_from_Lreuteri_530 = get_Lreu_from_tp(Lreu_seq,Lreuteri_530_seq,Lreuteri_530)

# %% <special step: process Lreu_from_iML1515 (have periplasm(_p) compartment)>

Lreu_from_iML1515 = My_def.model_refine.remove_compartment(Lreu_from_iML1515,compartment = '_p')


# %% <general step: add reaction 'from' notes and save models>

Lreu_from_iNF517.id = 'Lreu_from_iNF517'
Lreu_from_iBT721.id = 'Lreu_from_iBT721'
Lreu_from_iML1515.id = 'Lreu_from_iML1515'
Lreu_from_Lreuteri_530.id = 'Lreu_from_iML1515'

My_def.merge_model.note_model_from(Lreu_from_iNF517, ['iNF517','BBH'])
My_def.merge_model.note_model_from(Lreu_from_iBT721, ['iBT721','BBH'])
My_def.merge_model.note_model_from(Lreu_from_iML1515, ['iML1515','BBH'])
My_def.merge_model.note_model_from(Lreu_from_Lreuteri_530, ['Lreu_from_Lreuteri_530','BBH'])

cobra.io.save_json_model(Lreu_from_iNF517,'Lreu_from_iNF517.json',sort=True)
cobra.io.save_json_model(Lreu_from_iBT721,'Lreu_from_iBT721.json',sort=True)
cobra.io.save_json_model(Lreu_from_iML1515,'Lreu_from_iML1515.json',sort=True)
cobra.io.save_json_model(Lreu_from_Lreuteri_530,'Lreu_from_Lreuteri_530.json',sort=True)


# %% <general step: select a main model, iNF517>

os.chdir('../')
print('----- Mergeing draft models  -----')
Lreu_draft_1 = Lreu_from_iNF517.copy()
cobra.io.save_json_model(Lreu_draft_1,'Lreu_draft_1.json',sort='True')


# %% <general step: merge models, add unique reactions from other templates, iML1515 and iBT721>

#   option 1
#   by cobra function
#filed
#Lreu = Lreu.merge(Lreu_from_iML1515,inplace=False)
#Lreu = Lreu.merge(Lreu_from_iBT721,inplace=False)

#   option 2
# by def function


#Lreu_draft_2, report_df = My_def.merge_model.merge_draftmodels(Lreu_draft_1,Lreu_from_Lreuteri_530)
print('\033[1;31;47m'+'Draft model compare with iML1515 report (same reaid , differetnt equation)')
Lreu_draft_2, report_df_from_iML1515 = My_def.merge_model.merge_draftmodels(Lreu_draft_1,Lreu_from_iML1515)
print('\033[1;31;47m'+'Draft model compare with iML1515 report (same reaid , differetnt equation')
Lreu_draft_2, report_df_from_1BT721 = My_def.merge_model.merge_draftmodels(Lreu_draft_2,Lreu_from_iBT721)
print('\033[1;31;47m'+'Draft model compare with iML1515 report (same reaid , differetnt equation')
Lreu_draft_2, report_df_from_1BT721 = My_def.merge_model.merge_draftmodels(Lreu_draft_2,Lreu_from_Lreuteri_530)

print('\033[0;34;48m')

# %% <general step: add exchange transport and biomass reactions >
    # code not general, depends on template model ex tran bio reaction ids

print('----- add exchange transport and biomass reactions  -----')
Lreu_draft_2 = add_inf517_etb_rea(Lreu_draft_2,iNF517)
cobra.io.save_json_model(Lreu_draft_2,'Lreu_draft_2.json',sort='True')


# %% <general step: check iNF517 and draft model FBA results >

print('----- FBA result  -----')
Lreu_draft_3 = Lreu_draft_2.copy()
iNF517.solver = 'cplex'
Lreu_draft_3.solver = 'glpk'
# Bug: solver status is 'infeasible'

iNF517.objective = "BIOMASS_LRE"
print('iNF517 Biomass:',iNF517.optimize())

Lreu_draft_3.objective = "BIOMASS_LRE"
print('Lreu_draft_3 Biomass:',Lreu_draft_3.optimize())


# %% <special step: reless release constraint avoide infeasible fab results (remove midum update rates )>

for rea in Lreu_draft_3.reactions:
    if 'EX' in rea.id:
        if rea.lower_bound<=0 and rea.upper_bound <= 0:
            rea.upper_bound = 0.0
        elif rea.lower_bound >= 0 and rea.upper_bound >= 0:
            rea.lower_bound = 0.0

for rea in iNF517.reactions:
    if 'EX' in rea.id:
        if rea.lower_bound<=0 and rea.upper_bound <= 0:
            rea.upper_bound = 0.0
        elif rea.lower_bound >= 0 and rea.upper_bound >= 0:
            rea.lower_bound = 0.0


# %% <general step: gap fill >

print('----- Gap filling  -----')
#note, if still failed , need change other way to gap fill, such as biomass partly gapfill, check template FVA even FBA results.
#solution_biomass = cobra.flux_analysis.gapfill(Lreu_draft_2, iNF517)
#filed
# note: need to set integer_threshold=1e-10
solution_biomass_f = cobra.flux_analysis.gapfilling.GapFiller(Lreu_draft_3, iNF517, demand_reactions=False, integer_threshold=1e-10)
solution_biomass = solution_biomass_f.fill(iterations=1)[0]
biomass_gaps_set = set([i.id for i in solution_biomass])
print('biomass gaps number:',len(biomass_gaps_set))

reaset = set([i.id for i in Lreu_draft_3.reactions])
metset = set([i.id for i in Lreu_draft_3.metabolites])

for i in biomass_gaps_set :
    rea = iNF517.reactions.get_by_id(i)
    rea.notes['from'] = ['iNF517','gap']
    #Bug !!!
    Lreu_draft_3.add_reaction(rea)

Lreu_draft_3.objective = "BIOMASS_LRE"
print('Lreu_draft_3 Biomass:',Lreu_draft_3.optimize())
cobra.io.save_json_model(Lreu_draft_3,'Lreu_draft_3.json',sort='True')


print('=====  Done =====')




