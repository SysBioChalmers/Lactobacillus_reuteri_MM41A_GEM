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
#%%
keep_move_dic = {
                "MCOATA":"MACPMT",
                "DALTAL_LRE":"DALTAL",
                "KAS15":"kaasIII",
                "ALATRS":"ALATRS_1",
                "ASPTRS":"ASPTRS_1",
                "CYSTRS":"CYSTRS_1",
                "HISTRS":"HISTRS_1",
                "LEUTRS":"LEUTRS_1",
                "METTRS":"METTRS_1",
                "PHETRS":"PHETRS_1",
                "ARGTRS":"ARGTRS_1",
                "SERTRS":"SERTRS_1",
                "GLYTRS":"GLYTRS_1",
                "ILETRS":"ILETRS_1",
                "LYSTRS":"LYSTRS_1",
                "THRTRS":"THRTRS_1",
                "TRPTRS":"TRPTRS_1",
                "TYRTRS":"TYRTRS_1",
                "VALTRS":"VALTRS_1",
                "ASNTRS":"ASNTRS_1",

                "EAR40x":"BTMAT1",
                "PRFGS":"PRFGS_1",
                "PPA2":"EPPP",

                "LEUt2r":"LEUt2",

                "ARBt2":"ARAB_Lt",
                "ARBt2r":"ARBt2",
                "SERt3":"SERt2r",
                "SERt2r":"SERt3",
                #"THRt2r":"THRt3",


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
Lreuteri_530 = cobra.io.load_json_model('../Step_02_DraftModels/Template/template_models/Lreuteri_530_standlized.json')
Lreuteri_530.reactions.get_by_id('DM_2ahbut_c').id =  'EX_2ahbut_c'

# open that reactions
for rea in check_bounds_set:
    try:
        Lreu_draft_3_refined.reactions.get_by_id(rea).bounds = Lreuteri_530.reactions.get_by_id(rea).bounds
    except:
        Lreu_draft_3_refined.reactions.get_by_id(rea).remove_from_model()


#%% <change biomass and gapfill>


Lreuteri_530.objective = "BIOMASS"
print('Lreuteri_530 Biomass:',Lreuteri_530.optimize())

Lreu_draft_3_refined.objective = "BIOMASS_LRE"
print('Lreu_draft_3_refined Biomass_LRE:',Lreu_draft_3_refined.optimize())


Lreu_draft_3_refined_copy = Lreu_draft_3_refined.copy()
for rea in Lreuteri_530.exchanges:
    try :
        Lreu_draft_3_refined.reactions.get_by_id(rea.id).bounds  = rea.bounds
    except:
        try:
            rea.notes['from'] = list(set(rea.notes['from']+ ['Lreuteri_530']))
        except:
            rea.notes['from'] = ['Lreuteri_530']
        Lreu_draft_3_refined.add_reaction(rea)
Lreu_draft_3_refined.add_reaction(Lreuteri_530.reactions.get_by_id('EX_2ahbut_c'))


Lreu_draft_3_refined.add_metabolites(Lreuteri_530.metabolites.get_by_id('ocdctrACP_c'))
Lreu_draft_3_refined.add_metabolites(Lreuteri_530.metabolites.get_by_id('dag_LRE_c'))
Lreu_draft_3_refined.add_metabolites(Lreuteri_530.metabolites.get_by_id('dgdag_LRE_c'))


for rea in Lreu_draft_3_refined.reactions:
    if '_LRE' in rea.id or '_LRE' in rea.reaction:
        try:
            rea.reaction = Lreuteri_530.reactions.get_by_id(rea.id).reaction
            rea.notes['from'] = list(set(rea.notes['from']+ ['Lreuteri_530']))
        except:
            rea.remove_from_model()

# %% remove
reaids = ['EX_na1_e','EX_btn_e','EX_k_e','EX_mg2_e','EX_lac__D_e','EX_12ppd__R_e','EX_3hpp_e','EX_raffin_e','EX_ocdcea_e']
for reaid in reaids:
    Lreu_draft_3_refined.reactions.get_by_id(reaid).remove_from_model()

for rea in Lreuteri_530.reactions:
    if '_LRE' in rea.id or '_LRE' in rea.reaction:
        try:
            Lreu_draft_3_refined.reactions.get_by_id(rea.id)
        except:
            Lreu_draft_3_refined.add_reaction(rea)

Lreu_draft_3_refined.reactions.get_by_id('ATPM').bounds = Lreuteri_530.reactions.get_by_id('ATPM').bounds
Lreu_draft_3_refined.reactions.get_by_id('NNDMBRT').reaction = Lreuteri_530.reactions.get_by_id('NNDMBRT').reaction
Lreu_draft_3_refined.reactions.get_by_id('MCMAT8').reaction = Lreuteri_530.reactions.get_by_id('MCMAT8').reaction
Lreu_draft_3_refined.reactions.get_by_id('HDMAT7').reaction = Lreuteri_530.reactions.get_by_id('HDMAT7').reaction

# %%
Lreu_draft_3_refined.objective = "BIOMASS"
print('Lreu_draft_3_refined Biomass:',Lreu_draft_3_refined.optimize())


# %% <gapfill>
# solution_part_biomass = cobra.flux_analysis.gapfill(Lreu_draft_3_refined, Lreuteri_530, demand_reactions=False)


# %% <gapfill partly>: ['PYDAMt', 'ALATA_Lr', 'AOBUTDs' ,'BTNt2i','FA161tr','FA182tr','FA183tr','AACPS183']


gap_biomass_mets = ['LIP_LRE_c','LTAtotal_LRE_c','adeadocbl_c','btn_c','pydx5p_c']

gap_partly_set = set(['PYDAMt', 'ALATA_Lr', 'AOBUTDs' ,'BTNt2i','FA161tr','FA182tr','FA183tr','AACPS183'])
# gapfill unstable!!! here are the gap if the function run well
try:
    for gap_biomass_met in gap_biomass_mets:

        reatem = cobra.Reaction('objec')
        Lreu_draft_3_refined.add_reactions([reatem])
        Lreuteri_530.add_reactions([reatem])

        rea_str = gap_biomass_met + ' --> '
        Lreu_draft_3_refined.reactions.get_by_id('objec').reaction = rea_str
        Lreuteri_530.reactions.get_by_id('objec').reaction = rea_str

        Lreuteri_530.objective = "objec"
        print('Lreuteri_530:', Lreuteri_530.optimize())

        Lreu_draft_3_refined.objective = "objec"
        print('Lreu_draft_3_refined:',Lreu_draft_3_refined.optimize())
        if Lreu_draft_3_refined.optimize().objective_value < 1e-10:
            solution_part_biomass = cobra.flux_analysis.gapfill(Lreu_draft_3_refined, Lreuteri_530, demand_reactions=False)
            gap_partly_set = gap_partly_set | set([i.id for i in solution_part_biomass[0]])
    try:
        Lreu_draft_3_refined.reactions.get_by_id('objec').remove_from_model()
        Lreuteri_530.reactions.get_by_id('objec').remove_from_model()
    except:
        pass
except:
    pass
try:
    Lreu_draft_3_refined.reactions.get_by_id('objec').remove_from_model()
    Lreuteri_530.reactions.get_by_id('objec').remove_from_model()
except:
    pass

# %% <add gap>   biomass_gaps_set = set([i.id for i in solution_part_biomass[0]])


print('biomass gaps number:',len(gap_partly_set))
for gap in gap_partly_set:
    #Lreu_draft_3_refined.add_reaction()
    rea = Lreuteri_530.reactions.get_by_id(gap)
    rea.notes['from'] = ['Lreuteri_530','gap']
    Lreu_draft_3_refined.add_reaction(rea)
Lreu_draft_3_refined.objective = "BIOMASS"
print('Lreu_draft_3_refined Biomass:',Lreu_draft_3_refined.optimize())



# #%%
# pfba_solution = cobra.flux_analysis.pfba(Lreuteri_530)
# need_fluxes = pfba_solution.fluxes[abs(pfba_solution.fluxes) > 0]
# #%%
# gap_list = []
# for need_id in need_fluxes.index:
#     try:
#         if Lreu_draft_3_refined.reactions.get_by_id(need_id).reaction != Lreuteri_530.reactions.get_by_id(need_id).reaction:
#             print('equaction')
#             print(Lreuteri_530.reactions.get_by_id(need_id),need_fluxes[need_id])
#             print(Lreu_draft_3_refined.reactions.get_by_id(need_id))
#             pass
#         if need_fluxes[need_id] > Lreu_draft_3_refined.reactions.get_by_id(need_id).upper_bound or need_fluxes[need_id] < Lreu_draft_3_refined.reactions.get_by_id(need_id).lower_bound:
#             print('bounds')
#             print(Lreuteri_530.reactions.get_by_id(need_id),need_fluxes[need_id],Lreuteri_530.reactions.get_by_id(need_id).bounds)
#             print(Lreu_draft_3_refined.reactions.get_by_id(need_id),need_fluxes[need_id],Lreu_draft_3_refined.reactions.get_by_id(need_id).bounds)
#             #pass
#
#     except:
#         #Lreu_draft_3_refined.add_reaction(Lreuteri_530.reactions.get_by_id(need_id))
#
#         gap_list.append(need_id)
# # %%
# for i in  gap_list:
#     print(Lreuteri_530.reactions.get_by_id(i))
#     Lreu_draft_3_refined.add_reaction(Lreuteri_530.reactions.get_by_id(i))
#
# Lreu_draft_3_refined.objective = "BIOMASS"
# print('Lreu_draft_3_refined Biomass:',Lreu_draft_3_refined.optimize())




# %% <general step: add _missing tag remove remove_useless_mets and gens >

for gen in Lreu_draft_3_refined.genes:

    if  not gen.id.startswith('MBLCLPDI_') and not gen.id.endswith('_missing'):
        print(gen.id)
        #gen.id = gen.id+'_missing'
        reas = Lreu_draft_3_refined.genes.get_by_id(gen.id).reactions
        for rea in reas:
            rea.gene_reaction_rule = re.sub(gen.id+'(?!_missing)',gen.id+'_missing',rea.gene_reaction_rule)
# dead end
metdeadset = set()
readeadset = set()
for met in Lreu_draft_3_refined.metabolites:
    if len(met.reactions) <= 1  :
        for i in met.reactions:
            readeadset.add(i.id)
            if '_LRE' not in i.reaction:
                i.remove_from_model()
        met.remove_from_model()
        metdeadset.add(met.id)

Lreu_draft_3_refined = My_def.model_refine.remove_useless_mets(Lreu_draft_3_refined)
Lreu_draft_3_refined = My_def.model_refine.remove_useless_genes(Lreu_draft_3_refined)





# %% <general step: add annotation from bigg>
# TODO: add annotation

#load data
# import framed
# from framed.io.sbml import load_cbmodel
# from carveme import config
#
# model = load_cbmodel('../bigg_database/universe_draft.xml', flavor=config.get('sbml', 'default_flavor'))
# framed.save_sbml_model(model, '../bigg_database/universe_draft_cobra.xml',flavor='cobra')

#bigg_draft = cobra.io.read_sbml_model('../bigg_database/universe_draft_cobra.xml')
#bigg_draft2 = cobra.io.read_sbml_model('../bigg_database/universe_draft.xml')


# bigg_draft = cobra.io.read_sbml_model('../bigg_database/universe_draft.xml')
# bigg_draft_met_set = set([i.id for i in bigg_draft.metabolites])

# for met in Lreu_draft_3_refined.metabolites:
#      if met.id in bigg_draft_met_set:
#          bigg_met = bigg_draft.metabolites.get_by_id(met.id)
#          if met.charge != bigg_met.charge:
#              print(met,met.charge,bigg_met,bigg_met.charge)


with open('../bigg_database/universal_model.pickle', 'rb') as f:
         bigg_model = pickle.load(f)

for met in bigg_model.metabolites:

    ann_dic = {}

    if len(met.annotation)==0:
            continue

    elif type(met.annotation[0])=='str':
        ann_dic[i[0]] = [i[1]]

    else:
        for i in met.annotation:

            if i[0] not in ann_dic.keys():
                ann_dic[i[0]] = [i[1]]
            else:
                ann_dic[i[0]] = ann_dic[i[0]]+ [i[1]]

    met.annotation = ann_dic

for rea in bigg_model.reactions:

    ann_dic = {}

    if len(rea.annotation)==0:
            continue

    elif type(rea.annotation[0])=='str':
        ann_dic[i[0]] = [i[1]]

    else:
        for i in rea.annotation:

            if i[0] not in ann_dic.keys():
                ann_dic[i[0]] = [i[1]]
            else:
                ann_dic[i[0]] = ann_dic[i[0]]+ [i[1]]
    rea.annotation = ann_dic


iML1515 = cobra.io.read_sbml_model('../Initial_data/template_models/iML1515.xml')

# %% <add annotation>
no_an_mets = set()
no_an_reas = set()
for met in Lreu_draft_3_refined.metabolites:
    if met.charge =='':     #add charge
        try:
            met.charge = iML1515.metabolites.get_by_id(met.id).charge
        except:
            try:
                met.charge = Lreuteri_530.metabolites.get_by_id(met.id).charge
            except:
                print(met,'nocharge!!!')
                no_an_mets.add(met.id)
    if met.formula == '' or met.formula == 'X':     #add charge
        try:
            met.formula = iML1515.metabolites.get_by_id(met.id).formula
        except:
            try:
                met.formula = Lreuteri_530.metabolites.get_by_id(met.id).formula
            except:
                print(met,'no formula!!!')
                no_an_mets.add(met.id)
    if len(met.annotation) <1:      #add annotation
        try:
            met.annotation = iML1515.metabolites.get_by_id(met.id).annotation
        except:
            try:
                met.annotation = Lreuteri_530.metabolites.get_by_id(met.id).annotation
            except:
                try:
                    met.annotation = bigg_model.metabolites.get_by_id(met.id).annotation
                except:

                    print(met,'no annotation!!!')
                    no_an_mets.add(met.id)

for rea in Lreu_draft_3_refined.reactions:

    if len(rea.annotation) <1:      #add annotation
        try:
            rea.annotation = iML1515.reactions.get_by_id(rea.id).annotation
        except:
            try:
                rea.annotation = Lreuteri_530.reactions.get_by_id(rea.id).annotation
            except:

                try:
                    rea.annotation = bigg_model.reactions.get_by_id(rea.id).annotation
                except:
                    print(rea,'no annotation!!!')
                    no_an_mets.add(rea.id)
    try:
        rea.notes['from']
    except:
        print(rea,'no notes from')
        rea.notes['from'] = ['Lreuteri_530']
try:
    Lreu_draft_3_refined.reactions.get_by_id('objec').remove_from_model()
except:
    pass

# TODO genes annotation






# %% <save >

cobra.io.write_sbml_model(Lreu_draft_3_refined,'Lreu_draft_3_refined_0901.xml')
cobra.io.save_json_model(Lreu_draft_3_refined,'Lreu_draft_3_refined_0901.json',sort='True')
My_def.io_file.model2txt(Lreu_draft_3_refined, 'Lreu_draft_3_refined_0901.txt')


# %%
# cobra.io.write_sbml_model(Lreu_draft_3_refined,'Lreu_draft_3_refined_0901.xml')
#
# comd = ' memote report snapshot --filename "Lreu_draft_3_refined_0901.html" Lreu_draft_3_refined_0901.xml'
# os.system(comd)
#
