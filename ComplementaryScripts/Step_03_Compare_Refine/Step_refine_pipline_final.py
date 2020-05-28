#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 2019-08-20

"""Step_refine_pipline_final.py
:description : basic check .mass balance,annotation,notes
:param : 
:returns: 
:rtype: 
"""

import os
import pickle
import re

import cobra
import pandas as pd
from Bio import SeqIO

import My_def

os.chdir('../../ComplementaryData/Step_03_Compare_Refine/')
# %% os.chdir('ComplementaryData/Step_03_Compare_Refine/')
print('----- loading data -----')
Lreu_draft_3_refined = cobra.io.load_json_model('Lreu_draft_3_refined_part03.json')

Lreuteri_530 = cobra.io.load_json_model('../Step_02_DraftModels/Template/template_models/Lreuteri_530_standlized.json')
iNF517 = cobra.io.load_json_model('../Step_02_DraftModels/Template/template_models/iNF517_standlized.json')
iML1515 = cobra.io.load_json_model('../Step_02_DraftModels/Template/template_models/iML1515_standlized.json')

# %% <general step: mass balance>

iML1515_mets = set([i.id for i in iML1515.metabolites])
iNF517_mets = set([i.id for i in iNF517.metabolites])
Lreuteri_530_mets = set([i.id for i in Lreuteri_530.metabolites])

iML1515_reas = set([i.id for i in iML1515.reactions])
iNF517_reas = set([i.id for i in iNF517.reactions])
Lreuteri_530_reas = set([i.id for i in Lreuteri_530.reactions])

iML1515_reas_balanced = set([i.id for i in iML1515.reactions if i.check_mass_balance() == {}])
iNF517_reas_balanced = set([i.id for i in iNF517.reactions if i.check_mass_balance() == {}])
Lreuteri_530_reas_balanced = set([i.id for i in Lreuteri_530.reactions if i.check_mass_balance() == {}])

Lreu_draft_3_refined_reas = set([i.id for i in Lreu_draft_3_refined.reactions])
Lreu_draft_3_refined_reas_balanced = set([i.id for i in Lreu_draft_3_refined.reactions if i.check_mass_balance() == {}])

print('unbalanced reactions number:\t', len(Lreu_draft_3_refined_reas - Lreu_draft_3_refined_reas_balanced))

# <metabolites>
for met_i in Lreu_draft_3_refined.metabolites:
    if met_i.id in iML1515_mets:
        met_temp = iML1515.metabolites.get_by_id(met_i.id)
    elif met_i.id in Lreuteri_530_mets:
        met_temp = Lreuteri_530.metabolites.get_by_id(met_i.id)
    elif met_i.id in iNF517_mets:
        met_temp = iNF517.metabolites.get_by_id(met_i.id)
    else:
        # print(met_i)
        continue

    if met_i.formula != met_temp.formula:
        # print(met_i.id, met_i.formula, met_temp.formula)
        if met_temp.formula not in ['X', 'R', '']:
            Lreu_draft_3_refined.metabolites.get_by_id(met_i.id).formula = met_temp.formula
            Lreu_draft_3_refined.metabolites.get_by_id(met_i.id).charge = met_temp.charge
    elif met_i.charge != met_temp.charge:
        # print(met_i.id, met_i.charge, met_temp.charge)
        if met_temp.formula not in ['X', 'R', '']:
            Lreu_draft_3_refined.metabolites.get_by_id(met_i.id).charge = met_temp.charge

Lreu_draft_3_refined_reas_balanced = set([i.id for i in Lreu_draft_3_refined.reactions if i.check_mass_balance() == {}])
print('unbalanced reactions number:\t', len(Lreu_draft_3_refined_reas - Lreu_draft_3_refined_reas_balanced))

#  <reactions>
for rea_i in (Lreu_draft_3_refined_reas - Lreu_draft_3_refined_reas_balanced) & Lreuteri_530_reas_balanced:
    if Lreu_draft_3_refined.reactions.get_by_id(rea_i).reaction != Lreuteri_530.reactions.get_by_id(rea_i).reaction:
        # print(Lreu_draft_3_refined.reactions.get_by_id(rea_i))
        # print(Lreuteri_530.reactions.get_by_id(rea_i))
        if '_LRE' not in Lreu_draft_3_refined.reactions.get_by_id(rea_i).reaction:
            pass
            # Lreu_draft_3_refined.reactions.get_by_id(rea_i).reaction = Lreuteri_530.reactions.get_by_id(rea_i).reaction
    for met_i in Lreu_draft_3_refined.reactions.get_by_id(rea_i).metabolites.keys():
        if met_temp.formula not in ['X', 'R', '']:
            Lreu_draft_3_refined.metabolites.get_by_id(met_i.id).formula = Lreuteri_530.metabolites.get_by_id(
                met_i.id).formula
            Lreu_draft_3_refined.metabolites.get_by_id(met_i.id).charge = Lreuteri_530.metabolites.get_by_id(
                met_i.id).charge

Lreu_draft_3_refined_reas_balanced = set([i.id for i in Lreu_draft_3_refined.reactions if i.check_mass_balance() == {}])
print('unbalanced reactions number:\t', len(Lreu_draft_3_refined_reas - Lreu_draft_3_refined_reas_balanced))

for rea_i in (Lreu_draft_3_refined_reas - Lreu_draft_3_refined_reas_balanced) & iML1515_reas_balanced:
    if Lreu_draft_3_refined.reactions.get_by_id(rea_i).reaction != iML1515.reactions.get_by_id(rea_i).reaction:
        # print(Lreu_draft_3_refined.reactions.get_by_id(rea_i))
        # print(iML1515.reactions.get_by_id(rea_i))
        # if rea_i not in  ['MPTS']:
        Lreu_draft_3_refined.reactions.get_by_id(rea_i).reaction = iML1515.reactions.get_by_id(rea_i).reaction
    for met_i in Lreu_draft_3_refined.reactions.get_by_id(rea_i).metabolites.keys():
        try:
            Lreu_draft_3_refined.metabolites.get_by_id(met_i.id).formula = iML1515.metabolites.get_by_id(
                met_i.id).formula
            Lreu_draft_3_refined.metabolites.get_by_id(met_i.id).charge = iML1515.metabolites.get_by_id(met_i.id).charge
        except:
            pass

Lreu_draft_3_refined_reas_balanced = set([i.id for i in Lreu_draft_3_refined.reactions if i.check_mass_balance() == {}])
print('unbalanced reactions number:\t', len(Lreu_draft_3_refined_reas - Lreu_draft_3_refined_reas_balanced))

#  <hand check>manual check

for met_i in Lreu_draft_3_refined.metabolites:  # replace ACP:'X' to 'C11H21N2O7PRS'
    if 'acp_c' in met_i.id:
        # try:
        #     met_i.id = met_i.id.replace('acp_c','ACP_c')
        # except:
        #     new_id = met_i.id.replace('acp_c','ACP_c')
        #     Lreu_draft_3_refined = My_def.merge_model.merge_reactionsid(Lreu_draft_3_refined, new_id, met_i.id)
        formula = Lreu_draft_3_refined.metabolites.get_by_id(met_i.id).formula
        Lreu_draft_3_refined.metabolites.get_by_id(met_i.id).formula = formula.replace('X', 'C11H21N2O7PRS')
    elif 'ACP_c' in met_i.id:
        formula = Lreu_draft_3_refined.metabolites.get_by_id(met_i.id).formula
        Lreu_draft_3_refined.metabolites.get_by_id(met_i.id).formula = formula.replace('X', 'C11H21N2O7PRS')

Lreu_draft_3_refined.reactions.get_by_id('FMETTRS').reaction = Lreuteri_530.reactions.get_by_id('FMETTRS').reaction
Lreu_draft_3_refined.reactions.get_by_id(
    'PROTRS').reaction = 'atp_c + pro__L_c + trnapro_c --> amp_c + ppi_c + protrna_c'

refine_dic = {
    "2h3mp_c": ['C6H11O3', -1],
    "2h3mp_e": ['C6H11O3', -1],
    "2h3mb_c": ['C5H9O3', -1],
    "dimp_c": ['C10H12N4O7P', -1],
    "fol_e": ['C19H17N7O6', -2],
    "MGD_c": ['C20H20N10O14P2S2Mo', -2],
    "ppoh_e": ['C3H8O', 0],
    "PreZ_c": ['C10H11N5O7P', -1],
    "RTAala_c": ['C219H433N27O218P27', 25],
    "PGlac2_c": ['C40H62N7O22', -1],
    "glutrnagln_c": ['C5H7NO3R', 0],

    "pa_LRE_c": ['C37.36H67.59O8P', -2],
    "RNA_LRE_c": ['C9.52H9.91N3.75O6.58P1', -1],
    "PROT_LRE_c": ['C4.91H8.91N1.42O1.5S0.37', 1],
    "CPS_LRE_c": ['C24H47O27P', -2],
    "DNA_LRE_c": ['C9.81H11.31N3.69O6P', -1],
    "LIP_LRE_c": ['C41.68H77.65N0.46O10.17P', -1],
    "LTAtotal_LRE_c": ['C162.46H312.39N15.3O136.3P20', -5],
}

for met_i, formula in refine_dic.items():
    Lreu_draft_3_refined.metabolites.get_by_id(met_i).formula = formula[0]
    Lreu_draft_3_refined.metabolites.get_by_id(met_i).charge = formula[1]

Lreu_draft_3_refined_reas_unbalanced = set(
    [i.id for i in Lreu_draft_3_refined.reactions if (i.check_mass_balance() != {}) and ("EX" not in i.id)])
print('unbalanced reactions number:\t', len(Lreu_draft_3_refined_reas_unbalanced))
# removed_set = set(['PUACGAMS','HDMAT7','MCMAT8','GLUTRS2','ACKVAL','PTAVAL','ACKILE','PTAILE','SBTD','SBTD_L_1','NTPP10','ACKLEU','PTALEU','PYNP4'])
#
# for rea_i in list(Lreu_draft_3_refined_reas_unbalanced - removed_set):
#
#     rea = Lreu_draft_3_refined.reactions.get_by_id(rea_i)
#     # if '_LRE' in rea.reaction:
#     #     continue
#     if rea.check_mass_balance() !={}:
#         print('\n', rea)
#         print(rea.check_mass_balance())


Lreu_draft_3_refined.objective = "BIOMASS"
Lreu_draft_3_refined.optimize()
print('Lreu_draft_3_refined Biomass_LRE:', Lreu_draft_3_refined.optimize())

# %% <general step: metabolites and reactions annotation>

bigg = pickle.load(open('../bigg_database/universal_model.pickle', 'rb'))
bigg_mets = set([i.id for i in bigg.metabolites])
bigg_reas = set([i.id for i in bigg.reactions])
# Lreu_draft_3_refined = cobra.io.load_json_model('Lreu_draft_3_refined_part03.json')

for met_i in Lreu_draft_3_refined.metabolites:
    annotation_i = met_i.annotation
    if met_i.id in iML1515_mets:
        annotation_i = My_def.model_refine.merge_annotation(annotation_i,
                                                            iML1515.metabolites.get_by_id(met_i.id).annotation)
    if met_i.id in Lreuteri_530_mets:
        annotation_i = My_def.model_refine.merge_annotation(annotation_i,
                                                            Lreuteri_530.metabolites.get_by_id(met_i.id).annotation)
    if met_i.id in iNF517_mets:
        annotation_i = My_def.model_refine.merge_annotation(annotation_i,
                                                            iNF517.metabolites.get_by_id(met_i.id).annotation)
    if met_i.id in bigg_mets:
        annotation_temp = bigg.metabolites.get_by_id(met_i.id).annotation
        annotation_temp = My_def.model_refine.convert_annotation(annotation_temp)
        annotation_i = My_def.model_refine.merge_annotation(annotation_i, annotation_temp)

    if 'bigg.metabolite' not in annotation_i.keys():
        annotation_temp = {'bigg.metabolite': met_i.id.split('_')[0]}
        annotation_i = My_def.model_refine.merge_annotation(annotation_i, annotation_temp)
    if 'kegg.compound' in annotation_i.keys():
        annotation_temp = annotation_i['kegg.compound']
        if type(annotation_temp) == list:
            annotation_i['kegg.compound'] = [i for i in annotation_temp if i.startswith('C')]
    if 'sbo' not in annotation_i:
        if '_LRE' in met_i.id:
            # pass
            annotation_temp = {'sbo': ['SBO:0000649']}  # Biomass
        else:
            annotation_temp = {'sbo': ['SBO:0000247']}  # Simple chemical
        annotation_i = My_def.model_refine.merge_annotation(annotation_i, annotation_temp)

    Lreu_draft_3_refined.metabolites.get_by_id(met_i.id).annotation = annotation_i

for rea_i in Lreu_draft_3_refined.reactions:
    annotation_i = rea_i.annotation
    if rea_i.id in iML1515_reas:
        annotation_temp = iML1515.reactions.get_by_id(rea_i.id).annotation
        annotation_i = My_def.model_refine.merge_annotation(annotation_i, annotation_temp)
    if rea_i.id in Lreuteri_530_reas:
        annotation_temp = Lreuteri_530.reactions.get_by_id(rea_i.id).annotation
        annotation_i = My_def.model_refine.merge_annotation(annotation_i, annotation_temp)
    if rea_i.id in iNF517_reas:
        annotation_temp = iNF517.reactions.get_by_id(rea_i.id).annotation
        annotation_i = My_def.model_refine.merge_annotation(annotation_i, annotation_temp)
    if rea_i.id in bigg_reas:
        annotation_temp = bigg.reactions.get_by_id(rea_i.id).annotation
        annotation_temp = My_def.model_refine.convert_annotation(annotation_temp)
        annotation_i = My_def.model_refine.merge_annotation(annotation_i, annotation_temp)
    if 'bigg.reaction' not in annotation_i.keys():
        annotation_i = My_def.model_refine.merge_annotation(annotation_i, {'bigg.reaction': rea_i.id})
    if 'rhea' in annotation_i.keys():
        annotation_temp = annotation_i['rhea']
        if type(annotation_temp) == list:
            annotation_i['rhea'] = [i.split('#')[0] for i in annotation_temp]
    # if 'sbo' not in annotation_i:
    if rea_i.id.startswith('EX'):
        annotation_temp = {'sbo': ['SBO:0000627']}  # Exchange reactions
    elif 'BIOMASS' in rea_i.id or '_LRE' in rea_i.reaction:  # :
        annotation_temp = {'sbo': ['SBO:0000629']}  # Biomass reaction
    elif '_c' in rea_i.reaction and '_e' in rea_i.reaction:
        annotation_temp = {'sbo': ['SBO:0000655']}  # Transport reaction
    else:
        annotation_temp = {'sbo': ['SBO:0000176']}
    annotation_i = My_def.model_refine.merge_annotation(annotation_i, annotation_temp)

    Lreu_draft_3_refined.reactions.get_by_id(rea_i.id).annotation = annotation_i

for rea_i in Lreu_draft_3_refined.reactions:
    for i, j in rea_i.annotation.items():
        if type(j) == str:
            rea_i.annotation[i] = [j]


# get metacyc annotation important for subsystems.


def fill_specific_database(model_1, to_database):
    model = model_1.copy()
    id_rea_list = [rea_i.id for rea_i in model.reactions]

    metanetx_rea_list = []
    bigg_rea_list = []
    kegg_rea_list = []
    biocyc_rea_list = []
    for rea_i in model.reactions:
        if type(rea_i.annotation['bigg.reaction']) == str:
            rea_i.annotation['bigg.reaction'] = [rea_i.annotation['bigg.reaction']]
        else:
            bigg_rea_list = bigg_rea_list + rea_i.annotation['bigg.reaction']
        if 'metanetx.reaction' in rea_i.annotation.keys():
            metanetx_rea_list = metanetx_rea_list + rea_i.annotation['metanetx.reaction']
        if 'kegg.reaction' in rea_i.annotation.keys():
            kegg_rea_list = kegg_rea_list + rea_i.annotation['kegg.reaction']
        if 'biocyc' in rea_i.annotation.keys():
            biocyc_rea_list = biocyc_rea_list + rea_i.annotation['biocyc']

    biocyc_rea_list = [i.replace('META:', '') for i in biocyc_rea_list]
    targetlist_from_id, MNX_IDlist = My_def.mapIDsViaMNXref('rxns', id_rea_list, 'bigg', to_database)
    targetlist_from_meatnetx, MNX_IDlist = My_def.mapIDsViaMNXref('rxns', metanetx_rea_list, 'metanetx', to_database)
    targetlist_from_kegg, MNX_IDlist = My_def.mapIDsViaMNXref('rxns', kegg_rea_list, 'kegg', to_database)
    targetlist_from_biocyc, MNX_IDlist = My_def.mapIDsViaMNXref('rxns', biocyc_rea_list, 'metacyc', to_database)
    targetlist_from_bigg, MNX_IDlist = My_def.mapIDsViaMNXref('rxns', bigg_rea_list, 'bigg', to_database)

    if to_database == 'metacyc':
        annotation_key = 'biocyc'
    if to_database == 'kegg':
        annotation_key = 'kegg.reaction'

    for rea_i in model.reactions:
        if annotation_key in rea_i.annotation.keys():
            to_temp = rea_i.annotation[annotation_key]
            if rea_i.annotation[annotation_key] != ['']:
                continue

        to_temp = [targetlist_from_id[id_rea_list.index(rea_i.id)]]
        list_ann = ['metanetx.reaction', 'kegg.reaction', 'biocyc', 'bigg.reaction']
        list_ann.remove(annotation_key)
        for ann_i in list_ann:
            if ann_i in rea_i.annotation.keys():
                for i in rea_i.annotation[ann_i]:
                    if ann_i == 'metanetx.reaction':
                        to_temp_i = targetlist_from_meatnetx[metanetx_rea_list.index(i)]
                    elif ann_i == 'kegg.reaction':
                        to_temp_i = targetlist_from_kegg[kegg_rea_list.index(i)]
                    elif ann_i == 'biocyc':
                        i = i.replace('META:', '')
                        to_temp_i = targetlist_from_biocyc[biocyc_rea_list.index(i)]
                    elif ann_i == 'bigg.reaction':
                        to_temp_i = targetlist_from_bigg[bigg_rea_list.index(i)]

                    if type(to_temp_i) == str:
                        to_temp_i = [to_temp_i]
                    to_temp.extend(to_temp_i)
        to_temp = set(to_temp) - {''}
        if len(to_temp) != 0:
            # print(to_temp)
            rea_i.annotation[annotation_key] = list(to_temp)
    return model


to_database = 'metacyc'
Lreu_draft_3_refined = fill_specific_database(Lreu_draft_3_refined, to_database)

to_database = 'kegg'
Lreu_draft_3_refined = fill_specific_database(Lreu_draft_3_refined, to_database)

# %% <general step: subsystem annotation>
metacyc_subsystem = {}
for i in open('../Initial_data/MetaCyc_subSystems.csv', 'r'):
    i = i.replace('\n', '')
    temp = i.split('\t')
    if temp[1] != '':
        metacyc_subsystem[temp[0]] = temp[1:]

kegg_subsystem = {}
for i in open('../Initial_data/KEGG_subSystems.csv', 'r'):
    i = i.replace('\n', '')
    temp = i.split('\t')
    if temp[1] != '':
        kegg_subsystem[temp[0]] = temp[1:]

sub_count = 0
for rea_i in Lreu_draft_3_refined.reactions:
    subsystem = []
    if 'biocyc' in rea_i.annotation.keys():
        for i in rea_i.annotation['biocyc']:
            i = i.replace('META:', '')
            if i in metacyc_subsystem.keys():
                subsystem.extend(metacyc_subsystem[i])
    elif 'kegg.reaction' in rea_i.annotation.keys():
        for i in rea_i.annotation['kegg.reaction']:
            if i in kegg_subsystem.keys():
                subsystem.extend(kegg_subsystem[i])
    subsystem = set(subsystem) - {''}
    if len(subsystem) > 0:
        rea_i.annotation['subsystem'] = list(set(subsystem))
        sub_count += 1

# %% <general step: notes>
for rea_i in Lreu_draft_3_refined.reactions:
    if rea_i.notes == []:
        rea_i.notes = {}
    if 'from' in rea_i.notes.keys():
        if 'Lreu_from_Lreuteri_530' in rea_i.notes['from']:
            Lreu_draft_3_refined.reactions.get_by_id(rea_i.id).notes['from'].remove('Lreu_from_Lreuteri_530')
        if 'transport' in rea_i.notes['from']:
            Lreu_draft_3_refined.reactions.get_by_id(rea_i.id).notes['from'].remove('transport')
            Lreu_draft_3_refined.reactions.get_by_id(rea_i.id).notes['from'].append('transport_reaction')
        if 'EX' in rea_i.notes['from']:
            Lreu_draft_3_refined.reactions.get_by_id(rea_i.id).notes['from'].remove('EX')
            Lreu_draft_3_refined.reactions.get_by_id(rea_i.id).notes['from'].append('exchange_reaction')
    else:
        if rea_i.id in iML1515_reas:
            notes_temp = {'from': ['iML1515', 'gap']}
        elif rea_i.id in Lreuteri_530_reas:
            notes_temp = {'from': ['Lreuteri_530', 'gap']}
        elif rea_i.id in iNF517_reas:
            notes_temp = {'from': ['iNF517', 'gap']}
        else:
            notes_temp = {'from': ['gap']}
        Lreu_draft_3_refined.reactions.get_by_id(rea_i.id).notes = notes_temp

# %% <general step: genes annotation and add _missing tag to missing genes >
print('blast ing')  # find missing genes

'''
os.system('mkdir blast_tmps/')
Lreu_draft_3_refined_seq = '../Step_02_DraftModels/Lreuteri_biogaia_v03_2.faa'
iNF517_seq = '../Step_02_DraftModels/Template/template_seqs/iNF517.faa'
iML1515_seq = '../Step_02_DraftModels/Template/template_seqs/iML1515.faa'
Lreuteri_530_seq = '../Step_02_DraftModels/Template/template_seqs/Lreuteri_530.faa'
iBT721_seq = '../Step_02_DraftModels/Template/template_seqs/iBT721.faa'

qseq_file = Lreu_draft_3_refined_seq
mk_db = 'diamond makedb --in '+ qseq_file +' -d blast_tmps/qseq_db \n'
os.system(mk_db)
print('diamond blasting...')
options = ' --top 10 --more-sensitive '

df = pd.DataFrame()
for sseq_file in [iNF517_seq,iML1515_seq,Lreuteri_530_seq,iBT721_seq]:
    diamond_blastp_cmd ='diamond blastp -d blast_tmps/qseq_db -q ' + sseq_file + options +' -o ' +\
                        'blast_tmps/blast_result_s_in_q.csv --evalue 1 --outfmt 6 qseqid sseqid evalue pident length bitscore ppos qcovhsp \n'

    os.system(diamond_blastp_cmd)
    names = ['qseqid', 'sseqid', 'evalue', 'pident', 'length', 'bitscore', 'ppos','qcovs']
    df1 = pd.read_csv('blast_tmps/blast_result_s_in_q.csv', sep = '\t', names = names)
    df = df.append(df1)
df.drop_duplicates(subset='qseqid', keep='first')
df.index = df['qseqid']
os.system('rm -R blast_tmps/')
df.to_csv('blast_result.csv',sep = '\t')
'''

df = pd.read_csv('blast_result.csv', sep='\t', header=0)
df.index = df['qseqid']
df_dict = df['sseqid'].to_dict()

replace_gene_list = []
gene_notes = {}
for gen_i in Lreu_draft_3_refined.genes:  # find missing genes and notes
    old_tag = gen_i.id.replace('_missing', '')
    if old_tag in df_dict.keys():
        reas = list(gen_i.reactions.copy())
        new_id = df_dict[old_tag]
        # gen_i.id = new_id
        for rea_i in reas:
            rea_i.gene_reaction_rule = re.sub('_missing', '', rea_i.gene_reaction_rule)
            rea_i.gene_reaction_rule = re.sub(old_tag, new_id, rea_i.gene_reaction_rule)
        gene_notes[new_id] = {'old_tag': old_tag, 'old_tag_missing': old_tag + '_missing'}
    elif not gen_i.id.startswith('MBLCLPDI_') and not gen_i.id.endswith('_missing'):
        replace_gene_list.append(gen_i.id)
for gen_i in gene_notes.keys():
    Lreu_draft_3_refined.genes.get_by_id(gen_i).notes = gene_notes[new_id]

for gen_i in replace_gene_list:  # add missing
    rea_list = list(Lreu_draft_3_refined.genes.get_by_id(gen_i).reactions.copy())
    for rea_i in rea_list:
        rea_i.gene_reaction_rule = re.sub(gen_i + '(?!_missing)', gen_i + '_missing', rea_i.gene_reaction_rule)

Lreu_draft_3_refined = My_def.model_refine.remove_useless_genes(Lreu_draft_3_refined)

for rea_i in Lreu_draft_3_refined.reactions:
    if '_missing' in rea_i.gene_reaction_rule:
        # print(rea_i.gene_reaction_rule)
        rea_i.notes['missing_gene_reaction_rule'] = rea_i.gene_reaction_rule
        replace_list = []
        for gene_i in rea_i.genes:
            if '_missing' in gene_i.id:
                replace_list.append(gene_i.id)
        for i in replace_list:
            rea_i.gene_reaction_rule = rea_i.gene_reaction_rule.replace(i, 'check_rea_notes_for_genes_missing')
        # print(rea_i.gene_reaction_rule)

Lreu_draft_3_refined = My_def.model_refine.remove_useless_genes(Lreu_draft_3_refined)

print('genes annotation')

Lreu_draft_3_refined_gbk = '../Step_01_Sequences_analysis/Lreuteri_biogaia_v03/Lreuteri_biogaia_v03.gbk'
gen_set = [i.id for i in Lreu_draft_3_refined.genes]

open(Lreu_draft_3_refined_gbk, 'r')

locus_ncbiprotein_dict = dict()
locus_uniprot_dict = dict()
locus_refseq_dict = dict()
for seq in SeqIO.parse(Lreu_draft_3_refined_gbk, "gb"):
    for seq_feature in seq.features:

        if seq_feature.type == "CDS":
            locus_tag = seq_feature.qualifiers['locus_tag'][0]
            if 'similar to AA sequence:Lreuteri_refseq_v02.faa:' in seq_feature.qualifiers['inference'][-1]:
                ncbiprotein = seq_feature.qualifiers['inference'][-1].split(':')[-1]
                locus_ncbiprotein_dict[locus_tag] = ncbiprotein
            elif 'UniProtKB' in seq_feature.qualifiers['inference'][-1]:
                uniprot = seq_feature.qualifiers['inference'][-1].split(':')[-1]
                locus_uniprot_dict[locus_tag] = uniprot
uniprot_df = pd.read_csv('../Initial_data/uniprot-Lactobacillus+reuteri.tab', sep='\t', header=0)
uniprot_df = uniprot_df.fillna('')

for gen_i in Lreu_draft_3_refined.genes:
    if gen_i.id in locus_ncbiprotein_dict.keys():
        ncbiprotein = locus_ncbiprotein_dict[gen_i.id]
        gen_i.annotation['ncbiprotein'] = ncbiprotein
        df_temp = uniprot_df[uniprot_df['Cross-reference (RefSeq)'].str.match(ncbiprotein)]
        if df_temp.shape[0] > 0:
            gen_i.annotation['uniprot'] = [i for i in df_temp['Entry']]
            gen_i.annotation['ccds'] = [i.replace(';', '') for i in df_temp['Cross-reference (CCDS)'] if i != '']
            gen_i.annotation['kegg.genes'] = [i.replace(';', '') for i in df_temp['Cross-reference (KEGG)'] if i != '']
    elif gen_i.id in locus_uniprot_dict.keys():
        uniprot = locus_uniprot_dict[gen_i.id]
        gen_i.annotation['uniprot'] = uniprot
        df_temp = uniprot_df[uniprot_df['Entry'].str.match(uniprot)]
        if df_temp.shape[0] > 0:
            gen_i.annotation['ncbiprotein'] = [i.replace(';', '') for i in df_temp['Cross-reference (RefSeq)'] if
                                               i != '']
            gen_i.annotation['ccds'] = [i.replace(';', '') for i in df_temp['Cross-reference (CCDS)'] if i != '']
            gen_i.annotation['kegg.genes'] = [i.replace(';', '') for i in df_temp['Cross-reference (KEGG)'] if i != '']
    gen_i.annotation['refseq'] = ['NZ_ACGX02000007.1']

# %% <check duplactions>

print('\033[1;31;47m')
# case: general check duplaction
cofacters_set_1 = {'h2o_c', 'h2o_e', 'h_c', 'h_e'}
check_df = My_def.model_refine.check_duplicate_rea(Lreu_draft_3_refined, cofacters_set_1, remove=False)
check_id = list(check_df['id'])
print('Duplicate reactions: ')
for i in list(check_id):
    rea = Lreu_draft_3_refined.reactions.get_by_id(i)
    print(rea, rea.bounds, rea.gene_reaction_rule, rea.notes)
print('\033[0;34;48m')

# Manual Check: get the list from above report

# keep_move_dic{keep: remove}

keep_move_dic = {
    'THRt2r': 'THRt3',
}

for k, v in keep_move_dic.items():
    Lreu_draft_3_refined = My_def.merge_model.merge_reactionsid(Lreu_draft_3_refined, k, v)

# %% <general step: refine the model check cofactors > TODO : process cofactors!!
print('----- Refine cofacters ----- TODO!')

cofacters_set_2 = {'h2o_c', 'h2o_e', 'h_c', 'h_e',

                   'atp_c', 'adp_c', 'amp_c''pi_c', 'ppi_c', 'pi_e', 'ppi_e',
                   'gtp_c', 'gdp_c', 'gmp_c', 'itp_c', 'idp_c', 'imp_c', 'xtp_c', 'xmp_c', 'ditp_c', 'dimp_c'
                   }

cofacters_set_3 = {'h2o_c', 'h2o_e', 'h_c', 'h_e',

                   'nad_c', 'nadh_c', 'nadp_c', 'nadph_c', 'mqn7_c', 'mql7_c', 'mql8_c', 'mqn8_c',
                   'nadphx__R_c', 'nadphx__S_c', 'nadhx__R_c', 'nadhx__S_c',
                   'q8_c', 'q8h2_c', 'fad_c', 'fadh2_c',
                   'trdox_c', 'trdrd_c', 'grxox_c', 'grxrd_c', 'gthox_c', 'gthrd_c', 'pqqox_c', 'pqqrd_c',
                   }

print('\033[1;31;47m')
check_df = My_def.model_refine.check_duplicate_rea(Lreu_draft_3_refined, cofacters_set_1, remove=False)
# print(check_df)

check_id = list(check_df['id'])
print('Duplicate reactions: ')
for i in list(check_id):
    rea = Lreu_draft_3_refined.reactions.get_by_id(i)
    print(rea, rea.bounds, rea.gene_reaction_rule, rea.notes)
print('\033[0;34;48m')

# Manual Check: get the list from above report
# keep_move_dic{keep: remove}
keep_move_dic = {}

for k, v in keep_move_dic.items():
    Lreu_draft_3_refined = My_def.merge_model.merge_reactionsid(Lreu_draft_3_refined, k, v)

# %% <general step: check bounds == (0.0,0.0) >
check_bounds_set = set()
for i in Lreu_draft_3_refined.reactions:
    if i.bounds == (0.0, 0.0):
        check_bounds_set.add(i.id)
for rea_i in check_bounds_set:
    if rea_i in iML1515_reas:
        Lreu_draft_3_refined.reactions.get_by_id(rea_i).bounds = iML1515.reactions.get_by_id(rea_i).bounds

# <dead end process>


met_deadend_set = set()
rea_deadend_set = set()
gene_deadend_set = []
for met_i in Lreu_draft_3_refined.metabolites:
    if len(met_i.reactions) <= 1:
        met_deadend_set.add(met_i.id)
        for rea_i in met_i.reactions:
            rea_deadend_set.add(rea_i.id)
            for gen_i in rea_i.genes:
                gene_deadend_set.append(gen_i.id)

# remove:
removed = set([])
for rea_i in rea_deadend_set:
    remove_i = True
    for gen_i in Lreu_draft_3_refined.reactions.get_by_id(rea_i).genes:
        if 'MBLCLPDI_' in gen_i.id:
            if len(gen_i.reactions) == 1:
                remove_i = False
            # elif set([i.id for i in gen_i.reactions]) - rea_deadend_set == {}:
            #     remove_i = False
    if remove_i:
        Lreu_draft_3_refined.reactions.get_by_id(rea_i).remove_from_model()
        removed.add(rea_i)

Lreu_draft_3_refined = My_def.model_refine.remove_useless_genes(Lreu_draft_3_refined)
Lreu_draft_3_refined = My_def.model_refine.remove_useless_mets(Lreu_draft_3_refined)

# %% <growth check>
Lreu_draft_3_refined.objective = "BIOMASS"
Lreu_draft_3_refined.optimize()
print('Lreu_draft_3_refined Biomass_LRE:', Lreu_draft_3_refined.optimize())

# %% <save files>

for i in Lreu_draft_3_refined.metabolites:
    if i.compartment not in ['c', 'e']:
        i.compartment = i.id.split('_')[-1]

Lreu_draft_3_refined.id = 'iHL622'
Lreu_draft_3_refined.name = 'Lactobacillus reuteri ATCC PTA 6475 /(str. MM4-1A)'
Lreu_draft_3_refined.annotation = {'taxonomy': '548485', 'bigg.model': '',}
Lreu_draft_3_refined.notes = {'Model From': ['iNF517', 'LbReuteri', 'iBT721','iML1515']}


cobra.io.write_sbml_model(Lreu_draft_3_refined, '../../ModelFiles/iHL622.xml')
cobra.io.save_json_model(Lreu_draft_3_refined, '../../ModelFiles/iHL622.json', sort='True')
My_def.io_file.model2txt(Lreu_draft_3_refined, '../../ModelFiles/iHL622.txt')
cobra.io.save_matlab_model(Lreu_draft_3_refined, '../../ModelFiles/iHL622.mat')
cobra.io.save_yaml_model(Lreu_draft_3_refined, '../../ModelFiles/iHL622.yml')


# cobra.io.write_sbml_model(Lreu_draft_3_refined, '../../ModelFiles/Lreu_refined.xml')
# cobra.io.save_json_model(Lreu_draft_3_refined, '../../ModelFiles/Lreu_refined.json', sort='True')
# My_def.io_file.model2txt(Lreu_draft_3_refined, '../../ModelFiles/Lreu_refined.txt')
# cobra.io.save_matlab_model(Lreu_draft_3_refined, '../../ModelFiles/Lreu_refined.mat')
# cobra.io.save_yaml_model(Lreu_draft_3_refined, '../../ModelFiles/Lreu_refined.yml')

# memote report snapshot --filename "ModelFiles/Lreu_refined_memote_report_local.html" ModelFiles/iHL622.xml
# memote report snapshot --filename "ModelFiles/Lreu_refined_memote_report_local.html" ModelFiles/Lreu_refined.xml

# comd = 'memote report snapshot --filename "../../ModelFiles/Lreu_refined_memote_report_local.html" ../../ModelFiles/Lreu_refined.xml'
# os.system(comd)  # unstable , use wed app: https://memote.io/
