#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 2019-08-20

"""Step_refine_pipline.py
:description : script
:param : 
:returns: 
:rtype: 
"""

import My_def
import cobra
import os
import re

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


print('----- loading data -----')
os.chdir('../../ComplementaryData/Step_02_DraftModels/Template/')

# %% <general step: refine the model check duplications >

print('----- Refine duplications -----')
# case: special '_1' or '_2' in reaid compare and keep only one


Lreu_draft_3 = cobra.io.load_json_model('Lreu_draft_3.json')
Lreu_draft_3_refined = Lreu_draft_3.copy()
Lreu_draft_3_refined = remove_sepcical_dup(Lreu_draft_3_refined)


# case: general check duplaction
print('\033[1;31;47m')
check_df = My_def.model_refine.check_duplicate_rea(Lreu_draft_3_refined, cofacters_set = set(),remove = False )
#print(check_df)

check_id = list(check_df['id'])
print('Duplicate reactions: ')
for i in list(check_id):
    rea = Lreu_draft_3_refined.reactions.get_by_id(i)
    print(rea,rea.bounds,rea.gene_reaction_rule,rea.notes)
print('\033[0;34;48m')
# Manual Check: get the list from above report
# keep_move_dic{keep: remove}

keep_move_dic = {
                # "PRAIS" : "PRAIS_1",
                # "HEX1":"GLUK",
                # "GNK":"GNKr",
                # "CTPS1":"CTPS1__1",
                # "GLUCYS":"GLUCYSL",
                # "ATPM":"NTP1",
                # "FOLD3":"DHPS3",
                # "MPL":"MLTP4",
                # "LACZ":"GALSZ",
                # "PHETA1":"ASPTA6",
                # "ADNUC":"PNS1",
                # "PPGPPDP":"G35DP",
                # "GNNUC":"PNS2",
                # "INSH":"PNS3",
                # "XTSNH":"PNS4",
                # "URIH":"PYRNS1",
                # "G3PD1ir":"G3PD1",
                # "TDPDRR":"TDPDRR_1",
                # "GLUt2r":"GLUt6",
                # "DHPPDA":"DHPPDA__1",
                # "GALt2":"GALt",
                # "MALTt2":"MALTt",
                # "THRt2r":"THRt3",
                # "TYRt2r":"TYRt6",
                # "LEUt2r":"LEUt6",
                # "LYSt2r":"LYSt6",
                # "P5CR":"P5CRr",
                # "PTA2":"PDUL",
                # "PROTRS":"PROTRS_1",
                # "YUMPS":"PSUDS",
                # 'PROTS_LRE':'PROTS_LRE_v2',
                # 'MALt2r':'MALt6',
                                        }

for k,v in keep_move_dic.items():
    rea1 = Lreu_draft_3_refined.reactions.get_by_id(k)
    rea2 = Lreu_draft_3_refined.reactions.get_by_id(v)
    rea1.gene_reaction_rule = My_def.merge_model.merge_gprule(rea1.gene_reaction_rule, rea2.gene_reaction_rule)
    rea1.notes['from'] = list(set(rea1.notes['from']+rea2.notes['from']))
    Lreu_draft_3_refined.reactions.get_by_id(v).remove_from_model()



# %% <general step: refine the model check cofacters >

print('----- Refine cofacters -----')

cofacters_set_1 = { 'h2o_c', 'h2o_e', 'h_c', 'h_e'}


cofacters_set_2 = { 'h2o_c', 'h2o_e', 'h_c', 'h_e',

                    'atp_c', 'adp_c','amp_c''pi_c', 'ppi_c','pi_e','ppi_e',
                    'gtp_c','gdp_c','gmp_c','itp_c','idp_c','imp_c','xtp_c','xmp_c','ditp_c','dimp_c'
                    }

cofacters_set_3 = { 'h2o_c', 'h2o_e', 'h_c', 'h_e',

                    'nad_c', 'nadh_c','nadp_c','nadph_c','mqn7_c','mql7_c','mql8_c', 'mqn8_c',
                    'nadphx__R_c', 'nadphx__S_c', 'nadhx__R_c', 'nadhx__S_c',
                    'q8_c', 'q8h2_c', 'fad_c', 'fadh2_c',
                    'trdox_c', 'trdrd_c', 'grxox_c', 'grxrd_c','gthox_c','gthrd_c','pqqox_c', 'pqqrd_c',
                    }

print('\033[1;31;47m')
check_df = My_def.model_refine.check_duplicate_rea(Lreu_draft_3_refined, cofacters_set_1, remove = False  )
#print(check_df)

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
                }

for k,v in keep_move_dic.items():
    try:
        rea1 = Lreu_draft_3_refined.reactions.get_by_id(k)
    except:
        Lreu_draft_3_refined.reactions.get_by_id(v).id = k
        print(k +' not in model, just replase rea_id')
        continue
    try:
        rea2 = Lreu_draft_3_refined.reactions.get_by_id(v)
    except:
        print(v +' not in model')
        continue
    rea1.gene_reaction_rule = My_def.merge_model.merge_gprule(rea1.gene_reaction_rule, rea2.gene_reaction_rule)
    rea1.notes['from'] = list(set(rea1.notes['from']+rea2.notes['from']))
    Lreu_draft_3_refined.reactions.get_by_id(v).remove_from_model()




# %% <general step: check bounds == (0.0,0.0) >
check_bounds_set = set()
for i in Lreu_draft_3_refined.reactions:
    if i.bounds == (0.0,0.0):
        check_bounds_set.add(i.id)

# %% <general step: add _missing tag remove remove_useless_mets and gens >

for gen in Lreu_draft_3_refined.genes:

    if  not gen.id.startswith('MBLCLPDI_') and not gen.id.endswith('_missing'):
        print(gen.id)
        #gen.id = gen.id+'_missing'
        reas = Lreu_draft_3_refined.genes.get_by_id(gen.id).reactions
        for rea in reas:
            rea.gene_reaction_rule = re.sub(gen.id+'(?!_missing)',gen.id+'_missing',rea.gene_reaction_rule)

Lreu_draft_3_refined = My_def.model_refine.remove_useless_mets(Lreu_draft_3_refined)
Lreu_draft_3_refined = My_def.model_refine.remove_useless_genes(Lreu_draft_3_refined)


# %% <save >
cobra.io.save_json_model(Lreu_draft_3_refined,'Lreu_draft_3_refined_19-08-19.json',sort='True')
My_def.io_file.model2txt(Lreu_draft_3_refined, 'Lreu_draft_3_refined.txt')
