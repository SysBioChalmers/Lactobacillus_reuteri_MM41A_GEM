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
import pickle




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
Lreu_draft_3 = cobra.io.load_json_model('../Step_02_DraftModels/Lreu_draft_3.json')
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



# %% <general step: refine the model check cofacters >
# TODO : process cofacters!!

# print('----- Refine cofacters -----')
#
#
# cofacters_set_2 = { 'h2o_c', 'h2o_e', 'h_c', 'h_e',
#
#                     'atp_c', 'adp_c','amp_c''pi_c', 'ppi_c','pi_e','ppi_e',
#                     'gtp_c','gdp_c','gmp_c','itp_c','idp_c','imp_c','xtp_c','xmp_c','ditp_c','dimp_c'
#                     }
#
# cofacters_set_3 = { 'h2o_c', 'h2o_e', 'h_c', 'h_e',
#
#                     'nad_c', 'nadh_c','nadp_c','nadph_c','mqn7_c','mql7_c','mql8_c', 'mqn8_c',
#                     'nadphx__R_c', 'nadphx__S_c', 'nadhx__R_c', 'nadhx__S_c',
#                     'q8_c', 'q8h2_c', 'fad_c', 'fadh2_c',
#                     'trdox_c', 'trdrd_c', 'grxox_c', 'grxrd_c','gthox_c','gthrd_c','pqqox_c', 'pqqrd_c',
#                     }
#
# print('\033[1;31;47m')
# check_df = My_def.model_refine.check_duplicate_rea(Lreu_draft_3_refined, cofacters_set_1, remove = False  )
# #print(check_df)
#
# check_id = list(check_df['id'])
# print('Duplicate reactions: ')
# for i in list(check_id):
#     rea = Lreu_draft_3_refined.reactions.get_by_id(i)
#     print(rea,rea.bounds,rea.gene_reaction_rule,rea.notes)
# print('\033[0;34;48m')
#
# # Manual Check: get the list from above report
# # keep_move_dic{keep: remove}
#
# keep_move_dic = {
#                 # "LEUt2r":"LEUt2",
#                 #
#                 # "ARBt2":"ARAB_Lt",
#                 # "ARBt2r":"ARBt2",
#                 # "SERt3":"SERt2r",
#                 # "SERt2r":"SERt3",
#                 }
#
# for k,v in keep_move_dic.items():
#         Lreu_draft_3_refined = My_def.merge_model.merge_reactionsid(Lreu_draft_3_refined, k, v)
#



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

# %% <general step: add annotation from bigg>
# TODO: add annotation

# bigg_draft = cobra.io.read_sbml_model('../bigg_database/universe_draft.xml')
# bigg_draft_met_set = set([i.id for i in bigg_draft.metabolites])
#%%

# for met in Lreu_draft_3_refined.metabolites:
#      if met.id in bigg_draft_met_set:
#          bigg_met = bigg_draft.metabolites.get_by_id(met.id)
#          if met.charge != bigg_met.charge:
#              print(met,met.charge,bigg_met,bigg_met.charge)


# with open('../bigg_database/universal_model.pickle', 'rb') as f:
#         bigg_model = pickle.load(f)




#%% <biomass >
# realist = ['BIOMASS_LRE','BIOMASS_LRE_noATPnoH','DNAS_LRE','GALTAL',
#             'GAT1_LRE','LTAS_LRE','PAP_LRE','PROTS_LRE','PROTS_LRE_v3',
#             'RNAS_LRE','AGAT_LRE','CLPNS_LRE','DAGGT_LRE','DAGK_LRE',
#             'DASYN_LRE','PGSA_LRE','UGT1_LRE','UGT2_LRE','PGPP_LRE','LPGS_LRE',
#             'CPSS_LRE','DALTAL_LRE','LIP_LRE','GLTAL','DALTAL','LTA_total','BIOMASS']

Lreuteri_530 = cobra.io.load_json_model('../Step_02_DraftModels/Template/template_models/Lreuteri_530_standlized.json')
Lreuteri_530.reactions.get_by_id('DM_2ahbut_c').id =  'EX_2ahbut_c'

Lreuteri_530.objective = "BIOMASS"
print('Lreuteri_530 Biomass:',Lreuteri_530.optimize())


Lreu_draft_3_refined.objective = "BIOMASS_LRE"
print('Lreu_draft_3_refined Biomass_LRE:',Lreu_draft_3_refined.optimize())

# add new biomass
Lreu_draft_3_refined.add_reaction(Lreuteri_530.reactions.get_by_id('LIP_LRE'))
Lreu_draft_3_refined.add_reaction(Lreuteri_530.reactions.get_by_id('GLTAL'))
Lreu_draft_3_refined.add_reaction(Lreuteri_530.reactions.get_by_id('DALTAL'))
Lreu_draft_3_refined.add_reaction(Lreuteri_530.reactions.get_by_id('LTA_total'))
Lreu_draft_3_refined.add_reaction(Lreuteri_530.reactions.get_by_id('BIOMASS'))
Lreu_draft_3_refined.reactions.get_by_id('CPSS_LPL2').remove_from_model()


# add new transport
Lreu_draft_3_refined.add_reaction(Lreuteri_530.reactions.get_by_id('PIt2r'))
Lreu_draft_3_refined.add_reaction(Lreuteri_530.reactions.get_by_id('BTNt2i'))
Lreu_draft_3_refined.add_reaction(Lreuteri_530.reactions.get_by_id('D_LACt2'))
Lreu_draft_3_refined.add_reaction(Lreuteri_530.reactions.get_by_id('PYDAMt'))
Lreu_draft_3_refined.add_reaction(Lreuteri_530.reactions.get_by_id('FA161tr'))
Lreu_draft_3_refined.add_reaction(Lreuteri_530.reactions.get_by_id('FA181tr'))
Lreu_draft_3_refined.add_reaction(Lreuteri_530.reactions.get_by_id('FA182tr'))
Lreu_draft_3_refined.add_reaction(Lreuteri_530.reactions.get_by_id('FA183tr'))

Lreu_draft_3_refined.add_reaction(Lreuteri_530.reactions.get_by_id('DNMPPA'))
Lreu_draft_3_refined.add_reaction(Lreuteri_530.reactions.get_by_id('NH3c'))
Lreu_draft_3_refined.add_reaction(Lreuteri_530.reactions.get_by_id('PYDAMt'))
Lreu_draft_3_refined.add_reaction(Lreuteri_530.reactions.get_by_id('AACPS183'))
Lreu_draft_3_refined.add_reaction(Lreuteri_530.reactions.get_by_id('AOBUTDs'))
Lreu_draft_3_refined.add_reaction(Lreuteri_530.reactions.get_by_id('ALATA_Lr'))

for rea in Lreuteri_530.exchanges:
    try :
        Lreu_draft_3_refined.reactions.get_by_id(rea.id).bounds  = rea.bounds
    except:
        Lreu_draft_3_refined.add_reaction(rea)
Lreu_draft_3_refined.add_reaction(Lreuteri_530.reactions.get_by_id('EX_2ahbut_c'))

# refine
Lreu_draft_3_refined.reactions.get_by_id('THRt2r').reaction = Lreuteri_530.reactions.get_by_id('THRt2').reaction
Lreu_draft_3_refined.reactions.get_by_id('THRt2r').bounds = (-1000.0, 1000.0 )
Lreu_draft_3_refined.reactions.get_by_id('ASNt2r').bounds = (-1000.0, 1000.0 )
Lreu_draft_3_refined.reactions.get_by_id('ASPt2r').bounds = (-1000.0, 1000.0 )
Lreu_draft_3_refined.reactions.get_by_id('METt2r').bounds = (-1000.0, 1000.0 )
Lreu_draft_3_refined.reactions.get_by_id('PDH').bounds = (-1000.0, 1000.0 )
Lreu_draft_3_refined.reactions.get_by_id('NADH4').bounds = (-1000.0, 1000.0 )

Lreu_draft_3_refined.reactions.get_by_id('PROTS_LRE').reaction = Lreuteri_530.reactions.get_by_id('PROTS_LRE').reaction
Lreu_draft_3_refined.reactions.get_by_id('PROTS_LRE').bounds = Lreuteri_530.reactions.get_by_id('PROTS_LRE').bounds
# keep_move_dic = {
#                 "DALTAL_LPL":"DALTAL",
#                 "GLTAL_LPL":"GLTAL",
#                 }
#
# for k,v in keep_move_dic.items():
#     Lreu_draft_3_refined = My_def.merge_model.merge_reactionsid(Lreu_draft_3_refined, k, v)
# Lreu_draft_3_refined.add_reaction(Lreuteri_530.reactions.get_by_id('GLTAL'))
# Lreu_draft_3_refined.add_reaction(Lreuteri_530.reactions.get_by_id('DALTAL'))

Lreu_draft_3_refined.objective = "BIOMASS"
print('Lreu_draft_3_refined Biomass:',Lreu_draft_3_refined.optimize())

#Lreu_draft_3_refined_temp = Lreu_draft_3_refined.copy()



# for rea in Lreu_draft_3_refined.exchanges:
#     if rea.id.startswith('EX_'):
#         rea.bounds = (0.0, 0.0)



# Lreu_draft_3_refined.add_reaction(Lreuteri_530.reactions.get_by_id('DM_2ahbut_c'))
# Lreuteri_530.reactions.get_by_id('DM_2ahbut_c').id =  'EX_2ahbut_c'
# Lreu_draft_3_refined.reactions.get_by_id('DM_2ahbut_c').id =  'EX_2ahbut_c'
#solution_biomass = cobra.flux_analysis.gapfill(Lreu_draft_3_refined, Lreuteri_530, demand_reactions=False)


# gapfill
# Lreuteri_530.solver = 'cplex'
# Lreu_draft_3_refined.solver = 'glpk'
# solution_biomass_f = cobra.flux_analysis.gapfilling.GapFiller(Lreu_draft_3_refined, Lreuteri_530, demand_reactions=False, integer_threshold=1e-10)
# solution_biomass = solution_biomass_f.fill(iterations=1)[0]
#
pfba_solution = cobra.flux_analysis.pfba(Lreuteri_530)
need_fluxes = pfba_solution.fluxes[abs(pfba_solution.fluxes) > 0]

gap_list = []
for need_id in need_fluxes.index:
    try:
        if Lreu_draft_3_refined.reactions.get_by_id(need_id).reaction != Lreuteri_530.reactions.get_by_id(need_id).reaction:
            #print(Lreuteri_530.reactions.get_by_id(need_id))
            pass
        if need_fluxes[need_id] > Lreu_draft_3_refined.reactions.get_by_id(need_id).upper_bound or need_fluxes[need_id] < Lreu_draft_3_refined.reactions.get_by_id(need_id).lower_bound:
            print(Lreuteri_530.reactions.get_by_id(need_id),need_fluxes[need_id],Lreuteri_530.reactions.get_by_id(need_id).bounds)
            #pass

    except:
        Lreu_draft_3_refined.add_reaction(Lreuteri_530.reactions.get_by_id(need_id))

        #gap_list.append(need_id)
for i in  gap_list:
    print(Lreuteri_530.reactions.get_by_id(i))
    Lreu_draft_3_refined.add_reaction(Lreuteri_530.reactions.get_by_id(i))

Lreu_draft_3_refined.objective = "BIOMASS"
print('Lreu_draft_3_refined Biomass:',Lreu_draft_3_refined.optimize())


# biomass_gaps_set = set([i.id for i in solution_biomass])
# print('biomass gaps number:',len(biomass_gaps_set))
#
# reaset = set([i.id for i in Lreu_draft_3.reactions])
# metset = set([i.id for i in Lreu_draft_3.metabolites])
#
# for i in biomass_gaps_set :
#     rea = iNF517.reactions.get_by_id(i)
#     rea.notes['from'] = ['iNF517','gap']
#     #Bug !!!
#     Lreu_draft_3.add_reaction(rea)


# %% <save >
#cobra.io.save_json_model(Lreu_draft_3_refined,'Lreu_draft_3_refined_19-08-19.json',sort='True')
#My_def.io_file.model2txt(Lreu_draft_3_refined, 'Lreu_draft_3_refined.txt')
