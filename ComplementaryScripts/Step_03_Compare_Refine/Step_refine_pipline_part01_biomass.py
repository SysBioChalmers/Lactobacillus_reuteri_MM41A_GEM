#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 2019-08-20

"""Step_refine_pipline_final.py
:description : script to make the model have basic functions and overview of the model.
:param :
:returns:
:rtype:
"""

import os
import pickle
import re

import cobra

import My_def


def remove_sepcical_dup(Lreu_draft_3_refined):
    # case: reactions id end with '_1' or '_2', compare and keep only one

    # use:Lreu_draft_3_refined = remove_sepcical_dup(Lreu_draft_3_refined)

    model = Lreu_draft_3_refined.copy()
    for rea in model.reactions:

        if rea.id.endswith('_1') or rea.id.endswith('_2'):

            try:
                rea2 = Lreu_draft_3_refined.reactions.get_by_id(re.sub('_.$', '', rea.id))
                if rea2.reaction == rea.reaction:
                    # keep no _1 or _2 (rea2) one
                    rea2.gene_reaction_rule = My_def.merge_model.merge_gprule(rea.gene_reaction_rule,
                                                                              rea2.gene_reaction_rule)
                    rea2.notes['from'] = rea2.notes['from'] + rea.notes['from']
                    rea2.notes['from'] = list(set(rea2.notes['from']))

                    rea.remove_from_model()

                else:
                    # keep have _1 or _2 (rea) one
                    # print(rea)
                    # print(rea2)

                    # manual checked:just keep only one (iMl1515 one)
                    # pass
                    # GLUTRS_2: atp_c + glu__L_c + trnaglu_c --> amp_c + glutrna_c + h_c + ppi_c
                    # GLUTRS: atp_c + glu__L_c + trnaglu_c --> amp_c + glutrna_c + ppi_c
                    # PRAIS_1: atp_c + fpram_c --> adp_c + air_c + h_c + pi_c
                    # PRAIS: atp_c + fpram_c --> adp_c + air_c + 2.0 h_c + pi_c
                    # TDPDRR_1: dtdprmn_c + nadp_c <=> dtdp4d6dm_c + h_c + nadph_c
                    # TDPDRR: dtdp4d6dm_c + h_c + nadph_c --> dtdprmn_c + nadp_c

                    rea.gene_reaction_rule = My_def.merge_model.merge_gprule(rea.gene_reaction_rule,
                                                                             rea2.gene_reaction_rule)
                    rea.notes['from'] = rea2.notes['from'] + rea.notes['from']
                    rea.notes['from'] = list(set(rea.notes['from']))

                    if rea.bounds == (0, 0):
                        rea.bounds = rea2.bounds

                    rea2.remove_from_model()

            except  KeyError:
                pass
    return model


os.chdir('../../ComplementaryData/Step_03_Compare_Refine/')
# os.chdir('ComplementaryData/Step_03_Compare_Refine/')
# %%
print('----- loading data -----')
Lreu_draft_3 = cobra.io.load_json_model('../Step_02_DraftModels/Lreu_draft_3.json')
Lreu_draft_3_refined = Lreu_draft_3.copy()
Lreu_draft_3_refined.id = 'Lreu_draft_3_refined'
Lreu_draft_3_refined.objective = "BIOMASS_LRE"
print('Lreu_draft_3_refined Biomass_LRE:', Lreu_draft_3_refined.optimize())

Lreuteri_530 = cobra.io.load_json_model('../Step_02_DraftModels/Template/template_models/Lreuteri_530_standlized.json')
Lreuteri_530.reactions.get_by_id('DM_2ahbut_c').id = 'EX_2ahbut_c'
Lreuteri_530_copy = Lreuteri_530.copy()

# % <general step: refine biomass and gap fill>
Lreuteri_530.objective = "BIOMASS"
print('Lreuteri_530 Biomass:', Lreuteri_530.optimize())

Lreu_draft_3_refined.objective = "BIOMASS_LRE"
print('Lreu_draft_3_refined Biomass_LRE:', Lreu_draft_3_refined.optimize())

# <media>
for rea in Lreuteri_530.exchanges:  # add Lrue 530 exchange reactions
    try:
        Lreu_draft_3_refined.reactions.get_by_id(rea.id).bounds = rea.bounds
    except:
        try:
            rea.notes['from'] = list(set(rea.notes['from'] + ['Lreuteri_530', 'exchange_reactions']))
        except:
            rea.notes['from'] = ['Lreuteri_530', 'exchange_reactions']
        Lreu_draft_3_refined.add_reaction(rea)
Lreu_draft_3_refined.add_reaction(Lreuteri_530.reactions.get_by_id('EX_2ahbut_c'))

Lreu_draft_3_refined.objective = "BIOMASS_LRE"
print('Lreu_draft_3_refined changed media Biomass_LRE:', Lreu_draft_3_refined.optimize())

# <biomass reactions>
for rea in Lreu_draft_3_refined.reactions:
    if '_LRE' in rea.id or '_LRE' in rea.reaction:
        try:
            if rea.reaction != Lreuteri_530.reactions.get_by_id(rea.id).reaction:
                # print(rea)
                # print(Lreuteri_530.reactions.get_by_id(rea.id))
                rea.reaction = Lreuteri_530.reactions.get_by_id(rea.id).reaction
                rea.notes['from'] = list(set(rea.notes['from'] + ['Lreuteri_530']))
        except:
            rea.remove_from_model()
for rea in Lreuteri_530.reactions:
    if '_LRE' in rea.id or '_LRE' in rea.reaction:
        try:
            Lreu_draft_3_refined.reactions.get_by_id(rea.id)
        except:
            Lreu_draft_3_refined.add_reaction(rea)

Lreu_draft_3_refined.objective = "BIOMASS"
print('Lreu_draft_3_refined Biomass:', Lreu_draft_3_refined.optimize())

# <DNA>
Lreu_draft_3_refined.reactions.get_by_id(
    'DNAS_LRE').reaction = '1.37 atp_c + 0.306 datp_c + 0.194 dctp_c + 0.194 dgtp_c + 0.306 dttp_c + 1.37 h2o_c --> ' \
                           'DNA_LRE_c + 1.37 adp_c + 1.37 h_c + 1.37 pi_c + ppi_c '
# BOFdat {'dATP': -0.9788862507477464, 'dTTP': -1.00798745022546, 'dCTP': -0.6745686590765152, 
# 'dGTP': -0.5923042839569839, 'ppi_c': 3.2537466440067053} {'A': 0.3056069180215816, 'T': 0.3056069180215816, 
# 'C': 0.19439308197841845, 'G': 0.19439308197841845} 
# <others> 
Lreu_draft_3_refined.reactions.get_by_id(
    'LIP_LRE').reaction = 'h2o_c + ctp_c + pa_LRE_c + glyc3p_c+ 0.23 lystrna_c --> LIP_LRE_c + 0.02 glyc_c + 0.23 trnalys_c + pi_c + cmp_c + ppi_c'
Lreu_draft_3_refined.reactions.get_by_id(
    'LIP_LRE').gene_reaction_rule = '(MBLCLPDI_01949 or MBLCLPDI_00053) and MBLCLPDI_01680 and MBLCLPDI_00188 and MBLCLPDI_00571 and MBLCLPDI_00776'

Lreu_draft_3_refined.reactions.get_by_id(
    'AGAT_LRE').reaction = '0.5 glyc3p_c + 0.0211 cpocdacp_c + 0.0258 hdeACP_c + 0.075 ocdacp_c + 0.0418 ocdctrACP_c + 0.3519 ocdcyaACP_c + 0.1591 octeACP_c + 0.2566 palmACP_c + 0.0687 tdeacp_c --> ACP_c + 0.5 pa_LRE_c'
Lreu_draft_3_refined.reactions.get_by_id(
    'LTA_total').reaction = 'pa_LRE_c + 3.2 udpg_c + 35.3 atp_c + 21.0 h2o_c + 20.0 ctp_c + 20.0 glyc3p_c+ 15.3 ala__D_c --> LTAtotal_LRE_c + 20.0 cmp_c + 20.0 ppi_c + 35.3 adp_c + 23.2 h_c + 36.3 pi_c + 3.2 udp_c'
Lreu_draft_3_refined.reactions.get_by_id('LTA_total').gene_reaction_rule = 'MBLCLPDI_00825 and MBLCLPDI_01474'
Lreu_draft_3_refined.reactions.get_by_id('LTA_total').id = 'LTA_LRE'
Lreuteri_530 = Lreuteri_530_copy.copy()

Lreu_draft_3_refined.objective = "BIOMASS"
print('Lreu_draft_3_refined Biomass:', Lreu_draft_3_refined.optimize())

remove_list = ['CLPNS_LRE', 'PGPP_LRE', 'LPGS_LRE', 'PGSA_LRE', 'DASYN_LRE', 'DAGK_LRE', 'LTAS_LRE', 'DAGGT_LRE',
               'PAP_LRE', 'GLTAL', 'GAT1_LRE', 'DALTAL', 'UGT2_LRE']  # ]
for rea_i in remove_list:
    try:
        Lreu_draft_3_refined.reactions.get_by_id(rea_i).remove_from_model()
    except:
        pass

Lreu_draft_3_refined.reactions.get_by_id('ATPM').bounds = Lreuteri_530.reactions.get_by_id('ATPM').bounds
Lreu_draft_3_refined.reactions.get_by_id('NNDMBRT').reaction = Lreuteri_530.reactions.get_by_id('NNDMBRT').reaction
# Lreu_draft_3_refined.reactions.get_by_id('MCMAT8').reaction = Lreuteri_530.reactions.get_by_id('MCMAT8').reaction
# Lreu_draft_3_refined.reactions.get_by_id('HDMAT7').reaction = Lreuteri_530.reactions.get_by_id('HDMAT7').reaction

print('gap filling')
# gap_biomass_mets = ['LIP_LRE_c','LTAtotal_LRE_c','adeadocbl_c','btn_c','pydx5p_c']
gap_biomass_mets = [i.id for i, j in Lreu_draft_3_refined.reactions.get_by_id('BIOMASS').metabolites.items() if j < 0]
gap_biomass_mets = set(gap_biomass_mets) - {'atp_c', 'coa_c', 'h2o_c', 'nad_c', }

# gap_partly_set = set(['PYDAMt', 'ALATA_Lr', 'AOBUTDs' ,'BTNt2i','FA161tr','FA182tr','FA183tr','AACPS183'])
gaps_set = set([])

# gap fill unstable!!! here are the gap if the function run well
try:
    model_template = Lreuteri_530.copy()
    model_gap = Lreu_draft_3_refined.copy()
    model_template.objective = "BIOMASS"
    # print('Lreuteri_530:', model_template.optimize())
    model_gap.objective = "BIOMASS"
    # print('Lreu_draft_3_refined:',model_gap.optimize())

    if model_gap.optimize().objective_value < 1e-10:
        solution_part_biomass = cobra.flux_analysis.gapfill(model_gap, model_template, demand_reactions=False)
        gaps_set = gaps_set | set([i.id for i in solution_part_biomass[0]])
except:
    try:
        for gap_biomass_met in gap_biomass_mets:
            model_template = Lreuteri_530.copy()
            model_gap = Lreu_draft_3_refined.copy()

            rea_temp = cobra.Reaction('object')
            model_template.add_reactions([rea_temp])
            model_gap.add_reactions([rea_temp])

            rea_str = gap_biomass_met + ' --> '
            model_template.reactions.get_by_id('object').reaction = rea_str
            model_gap.reactions.get_by_id('object').reaction = rea_str

            # print('biomass apart optimize: ',gap_biomass_met)
            model_template.objective = "object"
            # print('Lreuteri_530:', model_template.optimize())
            model_gap.objective = "object"
            # print('Lreu_draft_3_refined:',model_gap.optimize())

            if model_gap.optimize().objective_value < 1e-10:
                solution_part_biomass = cobra.flux_analysis.gapfill(model_gap, model_template, demand_reactions=False)
                gaps_set = gaps_set | set([i.id for i in solution_part_biomass[0]])
    except:
        pass

'''
pfba_solution = cobra.flux_analysis.pfba(Lreuteri_530)
need_fluxes = pfba_solution.fluxes[abs(pfba_solution.fluxes) > 0]
#%%
gap_list = []
for need_id in need_fluxes.index:
    try:
        if Lreu_draft_3_refined.reactions.get_by_id(need_id).reaction != Lreuteri_530.reactions.get_by_id(need_id).reaction:
            print('equaction')
            print(Lreuteri_530.reactions.get_by_id(need_id),need_fluxes[need_id])
            print(Lreu_draft_3_refined.reactions.get_by_id(need_id))
            pass
        if need_fluxes[need_id] > Lreu_draft_3_refined.reactions.get_by_id(need_id).upper_bound or need_fluxes[need_id] < Lreu_draft_3_refined.reactions.get_by_id(need_id).lower_bound:
            print('bounds')
            print(Lreuteri_530.reactions.get_by_id(need_id),need_fluxes[need_id],Lreuteri_530.reactions.get_by_id(need_id).bounds)
            print(Lreu_draft_3_refined.reactions.get_by_id(need_id),need_fluxes[need_id],Lreu_draft_3_refined.reactions.get_by_id(need_id).bounds)
            #pass

    except:
        #Lreu_draft_3_refined.add_reaction(Lreuteri_530.reactions.get_by_id(need_id))

        gap_list.append(need_id)
# %%
for i in  gap_list:
    print(Lreuteri_530.reactions.get_by_id(i))
    Lreu_draft_3_refined.add_reaction(Lreuteri_530.reactions.get_by_id(i))

Lreu_draft_3_refined.objective = "BIOMASS"
print('Lreu_draft_3_refined Biomass:',Lreu_draft_3_refined.optimize())
'''

print('biomass gaps number:', len(gaps_set), gaps_set)

#  <add gaps>   biomass_gaps_set = set([i.id for i in solution_part_biomass[0]])
for gap in gaps_set:
    # Lreu_draft_3_refined.add_reaction()
    rea = Lreuteri_530.reactions.get_by_id(gap)
    rea.notes['from'] = ['Lreuteri_530', 'gap']
    Lreu_draft_3_refined.add_reaction(rea)
Lreu_draft_3_refined.objective = "BIOMASS"
print('Lreu_draft_3_refined Biomass:', Lreu_draft_3_refined.optimize())

for rea in Lreu_draft_3_refined.reactions:
    if '_LRE' in rea.id or '_LRE' in rea.reaction:
        print(rea.id, rea.reaction, rea.gene_reaction_rule)
    if '_LPL' in rea.id or '_LPL' in rea.reaction:
        print(rea.id, rea.reaction, rea.gene_reaction_rule)

Lreu_draft_3_refined.objective = "BIOMASS"
Lreu_draft_3_refined.optimize()

# %% <media>
temp_model = Lreu_draft_3_refined.copy()
remove_rea_ids = ['EX_2aeppn_e', 'EX_2h3mb_e', '2H3MBt',
                  # 'EX_2hxic__L_e',
                  'EX_2mba_e', '2MBAt6', 'EX_2mbald_e',
                  'EX_2mpa_e', '2MPAt6',
                  'EX_3mba_e', '3MBAt6', 'EX_3mbal_e',
                  'EX_acgam_e',
                  'EX_acmana_e',
                  'EX_bzal_e', 'BZALt',
                  'EX_cit_e', 'EX_fe3_e', 'EX_galt_e', 'EX_spmd_e', 'EX_zn2_e',
                  'SO4t2', 'RIBt2', 'METt2r', 'HISt2r', 'ABUTt2r',
                  'EX_na1_e', 'EX_k_e', 'EX_mg2_e', 'EX_lac__D_e',
                  'EX_12ppd__R_e', 'EX_3hpp_e', 'EX_raffin_e', 'EX_ocdcea_e']
# for rea_id in remove_rea_ids:
#     temp_model.reactions.get_by_id(rea_id).remove_from_model()
#     print(temp_model.optimize())

for rea_id in remove_rea_ids:
    try:
        Lreu_draft_3_refined.reactions.get_by_id(rea_id).remove_from_model()
    except:
        pass

# for rea_i in Lreu_draft_3_refined.reactions:
#     if '_c' in rea_i.reaction and '_e' in rea_i.reaction:
#         print(rea_i,rea_i.gene_reaction_rule,rea_i.notes)

Lreu_draft_3_refined.objective = "BIOMASS"
Lreu_draft_3_refined.optimize()

# <false gaps>
remove_list = set([])
for rea_i in Lreu_draft_3_refined.reactions:
    try:
        if rea_i.notes['from'] == ['iNF517', 'gap']:
            temp_bound = rea_i.bounds
            print(rea_i)
            rea_i.bounds = (0, 0)
            solution = Lreu_draft_3_refined.optimize()
            if solution.objective_value > 0.001:
                remove_list.add(rea_i.id)
            print('biomass: \t', Lreu_draft_3_refined.optimize().objective_value)
            rea_i.bounds = temp_bound
    except:
        pass
remove_list = remove_list - {'HSDy', 'G5SADs'}

for rea in remove_list:
    Lreu_draft_3_refined.reactions.get_by_id(rea).remove_from_model()

Lreu_draft_3_refined.objective = "BIOMASS"
print('Lreu_draft_3_refined Biomass:', Lreu_draft_3_refined.optimize())

# %% <general step: check and remove duplications >

print('----- Refine duplications -----')
# case: special '_1' or '_2' in rea_id compare and keep only one

Lreu_draft_3_refined = remove_sepcical_dup(Lreu_draft_3_refined)

# case: general check duplication

cofacters_set_1 = {'h2o_c', 'h2o_e', 'h_c', 'h_e'}
print('\033[1;31;47m')
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
    "MCOATA": "MACPMT",
    "DALTAL_LRE": "DALTAL",
    "KAS15": "kaasIII",
    "ALATRS": "ALATRS_1",
    "ASPTRS": "ASPTRS_1",
    "CYSTRS": "CYSTRS_1",
    "HISTRS": "HISTRS_1",
    "LEUTRS": "LEUTRS_1",
    "METTRS": "METTRS_1",
    "PHETRS": "PHETRS_1",
    "ARGTRS": "ARGTRS_1",
    "SERTRS": "SERTRS_1",
    "GLYTRS": "GLYTRS_1",
    "ILETRS": "ILETRS_1",
    "LYSTRS": "LYSTRS_1",
    "THRTRS": "THRTRS_1",
    "TRPTRS": "TRPTRS_1",
    "TYRTRS": "TYRTRS_1",
    "VALTRS": "VALTRS_1",
    "ASNTRS": "ASNTRS_1",

    "EAR40x": "BTMAT1",
    "PRFGS": "PRFGS_1",
    "PPA2": "EPPP",

    "LEUt2r": "LEUt2",

    "ARBt2": "ARAB_Lt",
    "ARBt2r": "ARBt2",
    "SERt3": "SERt2r",
    "SERt2r": "SERt3",
    # "THRt2r":"THRt3",

    "AGAT_LRE": "AGAT_LPL",
    "CLPNS_LRE": "CLPNS_LPL",
    "CPSS_LRE": "CPSS_LPL",
    # "CPSS_LRE":"CPSS_LPL2",
    "DAGGT_LRE": "DAGGT_LPL",
    "DALTAL_LPL": "DALTAL",
    "GLTAL_LPL": "GLTAL",
    "DASYN_LRE": "LPGS_LPL",
    "PGSA_LRE": "PGSA_LPL"
}

for k, v in keep_move_dic.items():
    Lreu_draft_3_refined = My_def.merge_model.merge_reactionsid(Lreu_draft_3_refined, k, v)

Lreu_draft_3_refined.objective = "BIOMASS"
print('Lreuteri_530 Biomass:', Lreuteri_530.optimize())

# %% <general step: check cofactors >
# TODO : process cofactors!!

# print('----- Refine cofactors -----')
#
#
# cofactors_set_2 = { 'h2o_c', 'h2o_e', 'h_c', 'h_e',
#
#                     'atp_c', 'adp_c','amp_c''pi_c', 'ppi_c','pi_e','ppi_e',
#                     'gtp_c','gdp_c','gmp_c','itp_c','idp_c','imp_c','xtp_c','xmp_c','ditp_c','dimp_c'
#                     }
#
# cofactors_set_3 = { 'h2o_c', 'h2o_e', 'h_c', 'h_e',
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
    if i.bounds == (0.0, 0.0):
        check_bounds_set.add(i.id)

# open that reactions
for rea in check_bounds_set:
    try:
        Lreu_draft_3_refined.reactions.get_by_id(rea).bounds = Lreuteri_530.reactions.get_by_id(rea).bounds
    except:
        Lreu_draft_3_refined.reactions.get_by_id(rea).remove_from_model()

Lreu_draft_3_refined.objective = "BIOMASS"
print('Lreuteri_530 Biomass:', Lreuteri_530.optimize())

# %% <general step: add _missing tag remove remove_useless_mets and gens >

for gen in Lreu_draft_3_refined.genes:

    if not gen.id.startswith('MBLCLPDI_') and not gen.id.endswith('_missing'):
        # print(gen.id)
        # gen.id = gen.id+'_missing'
        reas = Lreu_draft_3_refined.genes.get_by_id(gen.id).reactions
        for rea in reas:
            rea.gene_reaction_rule = re.sub(gen.id + '(?!_missing)', gen.id + '_missing', rea.gene_reaction_rule)

# dead end
met_dead_set = set()
rea_dead_set = set()
for met in Lreu_draft_3_refined.metabolites:

    if len(met.reactions) <= 1:
        for i in met.reactions:
            rea_dead_set.add(i.id)
            if len(i.gene_reaction_rule) == 0:
                print('remove', i)
                i.remove_from_model()
        # met.remove_from_model()
        met_dead_set.add(met.id)

Lreu_draft_3_refined = My_def.model_refine.remove_useless_mets(Lreu_draft_3_refined)
Lreu_draft_3_refined = My_def.model_refine.remove_useless_genes(Lreu_draft_3_refined)

Lreu_draft_3_refined.objective = "BIOMASS"
print('Lreu_draft_3_refined Biomass:', Lreu_draft_3_refined.optimize())

# %% <general step: add annotation from bigg>
bigg = pickle.load(open('../bigg_database/universal_model.pickle', 'rb'))
bigg_mets = set([i.id for i in bigg.metabolites])
bigg_reas = set([i.id for i in bigg.reactions])
Lreuteri_530 = cobra.io.load_json_model('../Step_02_DraftModels/Template/template_models/Lreuteri_530_standlized.json')
iNF517 = cobra.io.load_json_model('../Step_02_DraftModels/Template/template_models/iNF517_standlized.json')
iML1515 = cobra.io.load_json_model('../Step_02_DraftModels/Template/template_models/iML1515_standlized.json')

iML1515_mets = set([i.id for i in iML1515.metabolites])
iNF517_mets = set([i.id for i in iNF517.metabolites])
Lreuteri_530_mets = set([i.id for i in Lreuteri_530.metabolites])

iML1515_reas = set([i.id for i in iML1515.reactions])
iNF517_reas = set([i.id for i in iNF517.reactions])
Lreuteri_530_reas = set([i.id for i in Lreuteri_530.reactions])

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

# %% <general step: notes>
for rea_i in Lreu_draft_3_refined.reactions:
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
        notes_temp = set([])
        if rea_i.id in iML1515_reas:
            notes_temp = {'from': ['iML1515', 'gap']}
        elif rea_i.id in Lreuteri_530_reas:
            notes_temp = {'from': ['Lreuteri_530', 'gap']}
        elif rea_i.id in iNF517_reas:
            notes_temp = {'from': ['iNF517', 'gap']}
        Lreu_draft_3_refined.reactions.get_by_id(rea_i.id).notes = notes_temp

# %% <save >
Lreu_draft_3_refined.objective = "BIOMASS"
print('Lreu_draft_3_refined Biomass:', Lreu_draft_3_refined.optimize())

for i in Lreu_draft_3_refined.metabolites:
    if i.compartment not in ['c', 'e']:
        i.compartment = i.id.split('_')[-1]

cobra.io.write_sbml_model(Lreu_draft_3_refined, 'Lreu_draft_3_refined_part01.xml')
cobra.io.save_json_model(Lreu_draft_3_refined, 'Lreu_draft_3_refined_part01.json', sort='True')
My_def.io_file.model2txt(Lreu_draft_3_refined, 'Lreu_draft_3_refined_part01.txt', sort='True')

comd = ' memote report snapshot --filename "Lreu_draft_3_refined_part01.html" Lreu_draft_3_refined_part01.xml'
os.system(comd)
# comd = ' memote report snapshot --filename "Lreuteri_530.html" Lreuteri_530.xml'
# os.system(comd)
