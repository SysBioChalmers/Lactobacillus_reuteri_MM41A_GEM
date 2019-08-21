#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 2019-08-21

"""Step_37strains_refine_pipeline.py
:description : script
:param : 
:returns: 
:rtype: 
"""


import My_def
import cobra
import os
import re
import pickle

from Bio import SeqIO



def remove_sepcical_dup(Lreu_draft_3_refined):
    #case '_1' or '_2' in reaid compare and keep only one

    # use:Lreu_draft_3_refined = remove_sepcical_dup(Lreu_draft_3_refined)

    model = Lreu_draft_3_refined.copy()
    for rea in model.reactions:

        if rea.id.endswith('_1') or rea.id.endswith('_2'):

            try:
                rea2 = Lreu_draft_3_refined.reactions.get_by_id(re.sub('_.$','',rea.id))
                if rea2.reaction == rea.reaction:
                    # keep no _1 or _2 (rea2) one
                    rea2.gene_reaction_rule = My_def.merge_model.merge_gprule(rea.gene_reaction_rule, rea2.gene_reaction_rule)
                    rea2.notes['from'] = rea2.notes['from']+rea.notes['from']
                    rea2.notes['from'] = list(set(rea2.notes['from']))

                    rea.remove_from_model()

                else:
                    # keep have _1 or _2 (rea) one
                    #print(rea)
                    #print(rea2)

                    # manual checked:just keep only one (iMl1515 one)
                    #pass
                    # GLUTRS_2: atp_c + glu__L_c + trnaglu_c --> amp_c + glutrna_c + h_c + ppi_c
                    # GLUTRS: atp_c + glu__L_c + trnaglu_c --> amp_c + glutrna_c + ppi_c
                    # PRAIS_1: atp_c + fpram_c --> adp_c + air_c + h_c + pi_c
                    # PRAIS: atp_c + fpram_c --> adp_c + air_c + 2.0 h_c + pi_c
                    # TDPDRR_1: dtdprmn_c + nadp_c <=> dtdp4d6dm_c + h_c + nadph_c
                    # TDPDRR: dtdp4d6dm_c + h_c + nadph_c --> dtdprmn_c + nadp_c

                    rea.gene_reaction_rule = My_def.merge_model.merge_gprule(rea.gene_reaction_rule, rea2.gene_reaction_rule)
                    rea.notes['from'] = rea2.notes['from']+rea.notes['from']
                    rea.notes['from'] = list(set(rea.notes['from']))

                    if rea.bounds==(0,0):
                        rea.bounds = rea2.bounds

                    rea2.remove_from_model()

            except  KeyError:
                pass
    return model

os.chdir('../../ComplementaryData/Step_03_Compare_Refine/')
print('----- loading data -----')

name_list = ['100_23','20_2','3c6','ATCC53608','CF48_3A',
             'DSM200016','I5007','IRT','JCM1112',
             'LTH2584','LTH5448','MM2_3','MM4_1A','SD2112',
             'TD1','TMW1_112','TMW1_656','lpuph','mlc3','LR1',
             'LR10','LR11','LR12','LR13','LR14',
             'LR17','LR18','LR19','LR2','LR3',
             'LR4','LR6','LR7','LR8','LR9'
            ]
reacount = []
metcount = []
genecount = []
gene_nomissing_count = []

for sp_name in name_list:

    Lreu_draft_3 = cobra.io.load_json_model('../Step_02_DraftModels/other_Lreuteri_moldebuilding/'+ sp_name +'Lreu_draft_3.json')
    Lreu_draft_3_refined = Lreu_draft_3.copy()

    # %% <general step: refine the model check duplications >

    print('----- Refine duplications -----')
    # case: special '_1' or '_2' in reaid compare and keep only one

    Lreu_draft_3_refined = remove_sepcical_dup(Lreu_draft_3_refined)


    # case: general check duplaction

    cofacters_set_1 = { 'h2o_c', 'h2o_e', 'h_c', 'h_e'}
    print('\033[1;31;47m')
    check_df = My_def.model_refine.check_duplicate_rea(Lreu_draft_3_refined, cofacters_set_1,remove = False )

    check_id = list(check_df['id'])
    print('Duplicate reactions: ')
    for i in list(check_id):
        rea = Lreu_draft_3_refined.reactions.get_by_id(i)
        print(rea,rea.bounds,rea.gene_reaction_rule,rea.notes)

    print('\033[0;34;48m')

    # Manual Check: get the list from above report
    # keep_move_dic{keep: remove}

    keep_move_dic = {
                    "LEUt2r":"LEUt2",

                    "ARBt2":"ARAB_Lt",
                    "ARBt2r":"ARBt2",
                    "SERt3":"SERt2r",
                    "SERt2r":"SERt3",

                    "AGAT_LRE":"AGAT_LPL",
                    "CLPNS_LRE":"CLPNS_LPL",
                    "CPSS_LRE":"CPSS_LPL",
                    #"CPSS_LRE":"CPSS_LPL2",
                    "DAGGT_LRE":"DAGGT_LPL",
                    "DALTAL_LPL":"DALTAL",
                    "GLTAL_LPL":"GLTAL",
                    "DASYN_LRE":"LPGS_LPL",
                    "PGSA_LRE":"PGSA_LPL"
                    }

    for k,v in keep_move_dic.items():
        Lreu_draft_3_refined = My_def.merge_model.merge_reactionsid(Lreu_draft_3_refined, k, v)


    # %% <general step: check bounds == (0.0,0.0) >
    check_bounds_set = set()
    for i in Lreu_draft_3_refined.reactions:
        if i.bounds == (0.0,0.0):
            check_bounds_set.add(i.id)

    # %% <general step: add _missing tag remove remove_useless_mets and gens >


    faa_file = '../Step_02_DraftModels/other_Lreuteri_seq/'+sp_name + '.faa'
    geneset = set()
    for seq in SeqIO.parse(faa_file, "fasta"):
        geneset.add(seq.id)

    for gen in Lreu_draft_3_refined.genes:


        if   gen.id not in geneset and not gen.id.endswith('_missing'):
            #print(gen.id)
            #gen.id = gen.id+'_missing'
            reas = Lreu_draft_3_refined.genes.get_by_id(gen.id).reactions
            for rea in reas:
                rea.gene_reaction_rule = re.sub(gen.id+'(?!_missing)',gen.id+'_missing',rea.gene_reaction_rule)

    Lreu_draft_3_refined = My_def.model_refine.remove_useless_mets(Lreu_draft_3_refined)
    Lreu_draft_3_refined = My_def.model_refine.remove_useless_genes(Lreu_draft_3_refined)

    locals()[sp_name+'_reaset'] = set([ i.id for i in Lreu_draft_3_refined.reactions ])
    locals()[sp_name+'_metset'] = set([ i.id for i in Lreu_draft_3_refined.metabolites ])
    locals()[sp_name+'_geneset'] = set([ i.id for i in Lreu_draft_3_refined.genes ])
    locals()[sp_name+'_gene_no_missingset'] = set([ i for i in locals()[sp_name+'_geneset'] if '_missing' not in i ])

    reacount.append(len(locals()[sp_name+'_reaset']))
    metcount.append(len(locals()[sp_name+'_metset']))
    genecount.append(len(locals()[sp_name+'_geneset']))
    gene_nomissing_count.append(len(locals()[sp_name+'_gene_no_missingset']))



    cobra.io.save_json_model(Lreu_draft_3_refined,sp_name+'Lreu_draft_3_refined_19-08-21.json',sort='True')
    #My_def.io_file.model2txt(Lreu_draft_3_refined, 'Lreu_draft_3_refined.txt')

core_rea = locals()[name_list[0]+'_reaset']
pan_rea = locals()[name_list[0]+'_reaset']
core_met = locals()[name_list[0]+'_metset']
pan_met = locals()[name_list[0]+'_metset']

for sp_name in name_list:
    core_rea = core_rea & locals()[sp_name+'_reaset']
    pan_rea = pan_rea | locals()[sp_name+'_reaset']
    
    core_met = core_met & locals()[sp_name+'_metset']
    pan_met = pan_met | locals()[sp_name+'_metset']


#%%

carveme_reacount = []
carveme_metcount = []
carveme_genecount = []

for sp_name in name_list:


    carveme_model = cobra.io.read_sbml_model('../Step_02_DraftModels/other_Lreuteri_moldebuilding/CarveMe/'+ sp_name +'.xml')
    locals()[sp_name+'_carveme_reaset'] = set([ i.id for i in carveme_model.reactions ])
    locals()[sp_name+'_carveme_metset'] = set([ i.id for i in carveme_model.metabolites ])
    locals()[sp_name+'_carveme_geneset'] = set([ i.id for i in carveme_model.genes ])

    carveme_reacount.append(len(locals()[sp_name+'_reaset']))
    carveme_metcount.append(len(locals()[sp_name+'_metset']))
    carveme_genecount.append(len(locals()[sp_name+'_geneset']))

core_carveme_rea = locals()[name_list[0]+'_carveme_reaset']
pan_carveme_rea = locals()[name_list[0]+'_carveme_reaset']
core_carveme_met = locals()[name_list[0]+'_carveme_metset']
pan_carveme_met = locals()[name_list[0]+'_carveme_metset']


for sp_name in name_list:
    core_carveme_rea = core_carveme_rea & locals()[sp_name+'_carveme_reaset']
    pan_carveme_rea = pan_carveme_rea | locals()[sp_name+'_carveme_reaset']
    
    core_carveme_met = core_carveme_met & locals()[sp_name+'_carveme_metset']
    pan_carveme_met = pan_carveme_met | locals()[sp_name+'_carveme_metset']