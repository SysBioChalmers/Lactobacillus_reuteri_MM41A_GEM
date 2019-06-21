#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 2019-06-08

"""Step_03_add_EX_tran_biomass.py
:description : script
:param : draft model
:returns: gap filed model
:rtype: 
"""

import cobra
import os
import pandas as pd
from cobra.flux_analysis import gapfill
import re
import My_def
from cobra.flux_analysis import flux_variability_analysis

os.chdir('../../ComplementaryData/Step_03_Compare_Refine/')

#%% load data
iNF517 = cobra.io.load_json_model('iNF517.json')
iML1515 = cobra.io.load_json_model('iML1515.json')
Lreu_merged = cobra.io.load_json_model('Lreu_merged.json')
Lreu_merged.id = 'Lreu_merged'


#%% exchange reactions and biomass

reaset = set([i.id for i in Lreu_merged.reactions])
metset = set([i.id for i in Lreu_merged.metabolites])

for rea in iNF517.reactions:
    if ('EX_' in rea.id):
        # case exchange reactions
        if rea.id not in reaset:
            rea.notes['from'] = ['iNF517','EX']
            Lreu_merged.add_reaction(rea)
            reaset.add(rea.id)
    elif ('_c' in rea.reaction) and ('_e' in rea.reaction):
        # case transport
        if rea.id not in reaset:
            rea.notes['from'] = ['iNF517','transport']
            Lreu_merged.add_reaction(rea)
            reaset.add(rea.id)
    elif ('_LRE' in rea.id) or ('LRE_c' in rea.reaction):
        # case biomass
        #print(rea)
        if rea.functional:
            if rea.id not in reaset:
                rea.notes['from'] = ['iNF517','biomass']
                Lreu_merged.add_reaction(rea)
                reaset.add(rea.id)
    elif rea.id in ['ATPM'] and rea.id not in reaset:
        #case ATPM
        rea.notes['from'] = ['iNF517''atp']
        Lreu_merged.add_reaction(rea)
        reaset.add(rea.id)

# %% change update reates . avoide infeasible fab result
for rea in Lreu_merged.reactions:
    if 'EX' in rea.id:
        if rea.lower_bound<=0 and rea.upper_bound <= 0:
            rea.upper_bound = 0.0
        elif rea.lower_bound >= 0 and rea.upper_bound >= 0:
            rea.lower_bound = 0.0

# for rea in iNF517.reactions:
#     if 'EX' in rea.id:
#         if rea.lower_bound<=0 and rea.upper_bound <= 0:
#             rea.upper_bound = 0.0
#         elif rea.lower_bound >= 0 and rea.upper_bound >= 0:
#             rea.lower_bound = 0.0

# %% gapfill
reaset2 = reaset.copy()
metset2 = metset.copy()
Lreu_merged2 = Lreu_merged.copy()
iNF517.solver = 'cplex'
Lreu_merged.solver = 'glpk'

# %% case check ATP first, find gaps for ATP
Lreu_merged.objective = "ATPM"
print('ATPM:',Lreu_merged.optimize())
# result  = 1.027 (no gaps)

iNF517.objective = "ATPM"
print('ATPM:',iNF517.optimize())

# gapfill ATP
#solution_atp = gapfill(Lreu_merged, iNF517, demand_reactions=False)
##solution_atp = cobra.flux_analysis.gapfilling.GapFiller(Lreu_merged, iNF517, demand_reactions=False, integer_threshold=1e-9)
##solution_atp = solution_atp.fill(iterations=1)

# for rea in solution_atp[0]:
#     print(rea)
#
#     if rea.id not in reaset:
#         rea.notes['from'] = ['iNF517','gap','atp']
#         Lreu_merged.add_reaction(rea)
#         reaset.add(rea.id)
#
# Lreu_merged.objective = "ATPM"
# Lreu_merged.optimize()

# %% case whole biomass gapfill

Lreu_merged.objective = "BIOMASS_LRE"
print('Biomass:',Lreu_merged.optimize())

iNF517.objective = "BIOMASS_LRE"
print('Biomass:',iNF517.optimize())
#filed
#solution_biomass = gapfill(Lreu_merged, iNF517, demand_reactions=True)

solution_biomass_f = cobra.flux_analysis.gapfilling.GapFiller(Lreu_merged, iNF517, demand_reactions=False, integer_threshold=1e-10)
solution_biomass = solution_biomass_f.fill(iterations=1)[0]
biomass_gaps_set = set([i.id for i in solution_biomass])
print('biomass gaps:',solution_biomass)

# for rea in solution_biomass[0]:
#     print(rea)
#     if rea.id not in reaset:
#         rea.notes['from'] = ['iNF517','gap','biomass']
#         Lreu_merged.add_reaction(rea)
#         reaset.add(rea.id)



#%% case gapfill partly

#check models
for i in Lreu_merged.reactions:
    try:
        rea = iNF517.reactions.get_by_id(i.id)
        if i.bounds != rea.bounds:
            print('bounds different')
            print('Lreu',i,i.bounds)
            print('517',rea,rea.bounds)
        elif i.reaction != rea.reaction:
            print('mets different')
            print('Lreu',i,i.bounds)
            print('517',rea,rea.bounds)
    except:
        pass


# %% <option 1 >  get biomass reactions
biomass_reas = []
for rea in iNF517.reactions:
    if 'LRE' in rea.id:
        biomass_reas.append(rea.id)
        print('biomass', rea)
# %%
biomass_part_gaps_set = set()
Lreu_merged = Lreu_merged.copy()
for reaid in biomass_reas:
    Lreu_merged.objective = reaid
    f1 = Lreu_merged.optimize()
    iNF517.objective = reaid
    f2 = iNF517.optimize()
    print("---- Run %s ----" %reaid)
    print('Lreu_merged %s ae obj flux is: %d ' % (reaid, f1.objective_value))
    print('iNF517 %s ae obj flux is: %f ' % (reaid, f2.objective_value))

    if f1.objective_value > 0:
        print('no gap')
        rea_set = set()
    elif f2.objective_value <0.000001:
        print('iNF517 Error!')
    else:
        try:
            #solution_production = gapfill(Lreu_merged, iNF517, demand_reactions=True);
            solution_production = cobra.flux_analysis.gapfilling.GapFiller(Lreu_merged, iNF517, demand_reactions=False, integer_threshold=1e-13)
            solution_production = solution_production.fill(iterations=1)
            rea_set = set(i.id for i in solution_production[0])
            biomass_part_gaps_set = biomass_part_gaps_set|rea_set
        except cobra.exceptions.Infeasible:
            print('gapfilling optimization failed (infeasible)')
        except RuntimeError:
            print('RuntimeError: failed to validate gapfilled model, try lowering the integer_threshold')

# %% <option 2 >  get biomass composition filed
# biomass_dic = dict()
# for k,v in iNF517.reactions.get_by_id('BIOMASS_LRE').metabolites.items():
#     biomass_dic[k.id] = v
#
# def producton_gap_reas(gap_model,all_model,procuction,type = 'met'):
#
#     objrea = cobra.Reaction('objrea')
#     gap_model.add_reaction(objrea);
#     gap_model.reactions.get_by_id('objrea').reaction = procuction + ' --> '
#     gap_model.objective = "objrea"
#
#     f = gap_model.optimize()
#     print(f)
#     if f.objective_value > 0:
#         print('no gap')
#         rea_set = set()
#     else:
#
#         objrea_all = cobra.Reaction('objrea')
#         all_model.add_reaction(objrea_all);
#         all_model.reactions.get_by_id('objrea').reaction = procuction + ' --> '
#         all_model.objective = "objrea"
#
#         try:
#             solution_production = gapfill(gap_model, all_model, demand_reactions=False,exchange_reactions = True);
#             rea_set = set(i.id for i in solution_production[0])
#         except cobra.exceptions.Infeasible:
#             print('gapfilling optimization failed (infeasible)')
#             #print('try fva')
#             # f = flux_variability_analysis(all_model)
#             # rea_set = set([i for i in nesf.index])
#             rea_set = set()
#
#         #all_model.reactions.get_by_id('objrea').remove_from_model()
#
#     #gap_model.reactions.get_by_id('objrea').remove_from_model()
#     return rea_set
#
# biomass_part_gaps_set_2 = set()
#
# for met in biomass_dic.keys():
#
#     print(met)
#     met_rea_set = producton_gap_reas(Lreu_merged,iNF517,met)
#     print(met_rea_set - reaset)
#     biomass_part_gaps_set_2 = biomass_part_gaps_set_2|met_rea_set - reaset
# Lreu_merged.reactions.get_by_id('objrea').remove_from_model()
# iNF517.reactions.get_by_id('objrea').remove_from_model()

# for id in biomass_part_gaps_set_2:
#     if 'filed' not in id:
#         rea = iNF517.reactions.get_by_id(id)
#         rea.notes['from'] = ['iNF517','gap']
#         if rea.id not in reaset:
#             Lreu_merged.add_reaction(rea)
#             reaset.add(rea.id)



#%%  fva result
Lreu_merged.objective = "BIOMASS_LRE"
f = flux_variability_analysis(iNF517)
need_fva  = f[(f['minimum']>0.0001) | (f['maximum']< -0.0001)]
fva_gaps_set = set(need_fva.index) - reaset

# for need_id in need_fva.index:
#     #if need_id not in myrealist:
#     if need_id not in reaset:
#         rea = iNF517.reactions.get_by_id(need_id)
#         rea.notes['from'] = ['iNF517','gap']
#         Lreu_merged.add_reaction(rea)
#         print(rea)
#         reaset.add(rea.id)
#
# Lreu_merged.objective = "BIOMASS_LRE"
# Lreu_merged.optimize()



# %%  fba result

iNF517.objective = "BIOMASS_LRE"
solution_fba = iNF517.optimize()

need_fba = solution_fba.fluxes[abs(solution_fba.fluxes) > 1e-10]
fba_gaps_set = set(need_fba.index) - reaset

# for need_id in need_fba.index:
#     #if need_id not in myrealist:
#     if need_id not in reaset:
#         rea = iNF517.reactions.get_by_id(need_id)
#         rea.notes['from'] = ['iNF517','gap']
#         Lreu_merged.add_reaction(rea)
#         print(rea)
#         reaset.add(rea.id)


# %% change biomass according to refmodel  (L. reuteri JCM1112 )

#Lreu_merged.reactions.get_by_id('AGAT_LRE').reaction = '0.12 2chdeacp_c + 0.005 2ctdeacp_c + 0.32 2cocdacp_c + 0.01 agly3p_LRE_c + 0.25 cpocdacp_c + 0.26 hdeacp_c + 0.02 ocdacp_c + 0.03 tdeacp_c --> acp_c + 0.01 pa_LRE_c'
Lreu_merged.reactions.get_by_id('CPSS_LRE').reaction = 'dtdprmn_c + 5.0 h2o_c + 2.0 udpg_c + 2.0 udpgal_c <=> CPS_LRE_c + dtdp_c + 6.0 h_c + 3.0 udp_c + ump_c'
Lreu_merged.reactions.get_by_id('DALTAL_LRE').reaction = '0.01 LTA_LRE_c + 25.0 ala__D_c + 25.0 atp_c --> 0.01 LTAala_LRE_c + 25.0 adp_c + 25.0 pi_c'
#Lreu_merged.reactions.get_by_id('GAT1_LRE').reaction = '0.12 2chdeacp_c + 0.005 2ctdeacp_c + 0.32 2cocdacp_c + 0.25 cpocdacp_c + glyc3p_c + 0.295 hdeacp_c + 0.01 ocdacp_c + 0.09 tdeacp_c --> acp_c + 0.01 agly3p_LRE_c'
Lreu_merged.reactions.get_by_id('LTAS_LRE').reaction = 	'0.01 d12dg_LRE_c + 0.25 pg_LRE_c --> 0.01 12dgr_LRE_c + 0.25 LTA_LRE_c'
Lreu_merged.reactions.get_by_id('PROTS_LRE_v3').reaction = '0.125 alatrna_c + 0.04 argtrna_c + 0.06 asntrna_c + 0.06 asptrna_c + 0.306 atp_c + 0.011 cystrna_c + 0.083 glntrna_c + ' \
                                                       '0.023 glutrna_c + 0.084 glytrna_c + 2.0 gtp_c + 2.307 h2o_c + 0.017 histrna_c + 0.043 iletrna_c + 0.078 leutrna_c + 0.066 lystrna_c + ' \
                                                       '0.022 mettrna_c + 0.034 phetrna_c + 0.04 protrna_c + 0.056 sertrna_c + 0.064 thrtrna_c + 0.006 trptrna_c + 0.028 tyrtrna_c + 0.06 valtrna_c ' \
                                                       '--> 0.001 PROT_LRE_v3_c + 0.306 adp_c + 2.0 gdp_c + 2.306 h_c + 2.306 pi_c + 0.125 trnaala_c + 0.04 trnaarg_c + 0.06 trnaasn_c + 0.06 trnaasp_c + ' \
                                                       '0.011 trnacys_c + 0.083 trnagln_c + 0.023 trnaglu_c + 0.084 trnagly_c + 0.017 trnahis_c + 0.043 trnaile_c + 0.078 trnaleu_c + 0.066 trnalys_c + ' \
                                                       '0.022 trnamet_c + 0.034 trnaphe_c + 0.04 trnapro_c + 0.056 trnaser_c + 0.064 trnathr_c + 0.006 trnatrp_c + 0.028 trnatyr_c + 0.06 trnaval_c'



# %% compare gaps add gap reactions:

print( len(biomass_gaps_set))
print( len(biomass_part_gaps_set))
#print( len(biomass_part_gaps_set_2))
print( len(fva_gaps_set))
print( len(fba_gaps_set))

# add reactions
for i in biomass_gaps_set:
    rea = iNF517.reactions.get_by_id(i)
    rea.notes['from'] = ['iNF517','gap','biomass']
    Lreu_merged.add_reaction(rea)
    #reaset.add(rea.id)

gset = set()
for i in Lreu_merged.genes:
    if 'MBLCL' in i.id:
        gset.add(i.id)
print(len(gset))


Lreu_merged.objective = "BIOMASS_LRE"
print(Lreu_merged.optimize())

# %% <save models>
Lreu_merged.objective = "BIOMASS_LRE"
cobra.io.save_json_model(Lreu_merged,'../Step_03_Compare_Refine/Lreu_merged_gapfiled.json')








