#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 2019-05-24

"""compare_bulid_model_from_template.py
:description : script
:param : 
:returns: 
:rtype: 
"""
import os
import cobra
import pandas as pd
import Step_03_build_model_from_template
import Step_03_template_modle
import Step_02_template_blast
import My_def


# Modeling pipeline

# load template model
#template_model = cobra.io.read_sbml_model('iJO1366.xml')
os.chdir('../../ComplementaryData/Step_02_DraftModels/Template/')
template_model = cobra.io.load_json_model('template_models/iML1515_standlized.json')
qseq_file = '../Lreuteri_biogaia_v03.faa'
sseq_file = 'template_seqs/' + 'iML1515' + '.faa'

# do blast
Step_03_build_model_from_template.do_blast(qseq_file, sseq_file)

BBHs = Step_03_build_model_from_template.extract_BBHs(1e-20,True,coverage=45)

candidate_rxns = Step_03_build_model_from_template.get_all_rxns_in_BBH(template_model, BBHs)

model1  = Step_03_build_model_from_template.build_model_from_template(candidate_rxns,BBHs,True)




blast_result_s_in_q, blast_result_q_in_s = Step_02_template_blast.blastp_pairwise(qseq_file, sseq_file, out_dir='')

#blast_result_s_in_q = 'blast/blast_result_s_in_q.csv'
#blast_result_q_in_s = 'blast/blast_result_q_in_s.csv'
blast_result_df = My_def.select_blast(blast_result_s_in_q, blast_result_q_in_s, best_match = True, evalue=1e-20,
                                      pident=0, length=0,
                                      bitscore=0, ppos=0, qcovs=45)#(bitscore=100, ppos=45)
model2 = Step_03_template_modle.get_model_from_template(template_model, blast_result_df,remove_missing_genes = False)


# blast_result_df2 = My_def.select_blast(blast_result_s_in_q, blast_result_q_in_s, best_match=False, evalue=1e-20,
#                                       pident=0, length=0,
#                                       bitscore=0, ppos=0, qcovs=45)#(bitscore=100, ppos=45)
# model3 = Step_03_template_modle.get_model_from_template(template_model, blast_result_df2)

# print (blast_result_df[blast_result_df['sseqid'] == 'MBLCLPDI_00092'])
# print (blast_result_df[blast_result_df['qseqid'] == 'lp_1729'])

geneset1 = set([i.id for i in model1.genes])
geneset2 = set([i.id for i in model2.genes])


reaset1 = set([i.id for i in model1.reactions])
reaset2 = set([i.id for i in model2.reactions])


BBHs2 = []
for a in blast_result_df[['sseqid','qseqid']].itertuples(index=False):
    a = tuple(a)
    BBHs2.append(a)
BBHs_set = set(BBHs)
BBHs2_set = set(BBHs2)
print('BBHs different: BBHs2_set - BBHs_set', BBHs2_set - BBHs_set)
print('BBHs different: BBHs_set - BBHs2_set', BBHs_set - BBHs2_set)


for i in reaset1- reaset2:
    print(model1.reactions.get_by_id(i).gene_reaction_rule)



candidate_rxns2 = Step_03_build_model_from_template.get_all_rxns_in_BBH(template_model, BBHs2)
model3  = Step_03_build_model_from_template.build_model_from_template(candidate_rxns2,BBHs2,True)


model4 = Step_03_template_modle.get_model_from_template(template_model, blast_result_df,remove_missing_genes = False)
geneset3 = set([i.id for i in model3.genes])
reaset3 = set([i.id for i in model3.reactions])

geneset4 = set([i.id for i in model4.genes])
reaset4 = set([i.id for i in model4.reactions])

reaset3 - reaset4


