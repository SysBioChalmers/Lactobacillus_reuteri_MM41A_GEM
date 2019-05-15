#!/usr/bin/env python
#-*- coding: utf-8 -*-
# Created by lhao at 2019-04-22
import os
import cobra
import re
import My_def
import pandas as pd
import pickle
from cobra.util.solver import linear_reaction_coefficients
from cobra import Reaction
from cobra.flux_analysis import gapfill
from cobra.flux_analysis import flux_variability_analysis

os.chdir('/Users/lhao/Box Sync/Project_Lreuteri/03_Compare_Refine/')

iNF517 = cobra.io.load_json_model('iNF517_initial.json')
Lreu_ca_gp = cobra.io.load_json_model('Lreu_ca_gp.json')
mymodel = cobra.io.load_json_model('add_two.json')

for rea in iNF517.reactions:
    if ('EX_' in rea.id):
        #exchange reactions
        rea.notes['from'] = ['iNF517']
        mymodel.add_reaction(rea)
    if ('_LLA' in rea.id) or ('LLA_c' in rea.reaction):
        if rea.functional:
            rea.notes['from'] = ['iNF517']
            mymodel.add_reaction(rea)

#mymodel.reactions.get_by_id('AGAT_LLA').reaction = '0.12 2chdeacp_c + 0.005 2ctdeacp_c + 0.32 2cocdacp_c + 0.01 agly3p_LLA_c + 0.25 cpocdacp_c + 0.26 hdeacp_c + 0.02 ocdacp_c + 0.03 tdeacp_c --> acp_c + 0.01 pa_LLA_c'
mymodel.reactions.get_by_id('CPSS_LLA').reaction = 'dtdp6dm_c + 5.0 h2o_c + 2.0 udpg_c + 2.0 udpgal_c <=> CPS_LLA_c + dtdp_c + 6.0 h_c + 3.0 udp_c + ump_c'
mymodel.reactions.get_by_id('DALTAL_LLA').reaction = '0.01 LTA_LLA_c + 25.0 ala__D_c + 25.0 atp_c --> 0.01 LTAala_LLA_c + 25.0 adp_c + 25.0 pi_c'
#mymodel.reactions.get_by_id('GAT1_LLA').reaction = '0.12 2chdeacp_c + 0.005 2ctdeacp_c + 0.32 2cocdacp_c + 0.25 cpocdacp_c + glyc3p_c + 0.295 hdeacp_c + 0.01 ocdacp_c + 0.09 tdeacp_c --> acp_c + 0.01 agly3p_LLA_c'
mymodel.reactions.get_by_id('LTAS_LLA').reaction = 	'0.01 d12dg_LLA_c + 0.25 pg_LLA_c --> 0.01 12dgr_LLA_c + 0.25 LTA_LLA_c'
mymodel.reactions.get_by_id('PROTS_LLA_v3').reaction = '0.125 alatrna_c + 0.04 argtrna_c + 0.06 asntrna_c + 0.06 asptrna_c + 0.306 atp_c + 0.011 cystrna_c + 0.083 glntrna_c + ' \
                                                       '0.023 glutrna_c + 0.084 glytrna_c + 2.0 gtp_c + 2.307 h2o_c + 0.017 histrna_c + 0.043 iletrna_c + 0.078 leutrna_c + 0.066 lystrna_c + ' \
                                                       '0.022 mettrna_c + 0.034 phetrna_c + 0.04 protrna_c + 0.056 sertrna_c + 0.064 thrtrna_c + 0.006 trptrna_c + 0.028 tyrtrna_c + 0.06 valtrna_c ' \
                                                       '--> 0.001 PROT_LLA_v3_c + 0.306 adp_c + 2.0 gdp_c + 2.306 h_c + 2.306 pi_c + 0.125 trnaala_c + 0.04 trnaarg_c + 0.06 trnaasn_c + 0.06 trnaasp_c + ' \
                                                       '0.011 trnacys_c + 0.083 trnagln_c + 0.023 trnaglu_c + 0.084 trnagly_c + 0.017 trnahis_c + 0.043 trnaile_c + 0.078 trnaleu_c + 0.066 trnalys_c + ' \
                                                       '0.022 trnamet_c + 0.034 trnaphe_c + 0.04 trnapro_c + 0.056 trnaser_c + 0.064 trnathr_c + 0.006 trnatrp_c + 0.028 trnatyr_c + 0.06 trnaval_c'

#cobra.io.write_sbml_model(mymodel,'mymodel.txt')
'''
for rea in mymodel.reactions:
    if '_LLA' in rea.id:
        rea.id = rea.id.replace('_LLA','_LRE')
        rea.name = rea.name.replace('_LLA','_LRE')
for met in mymodel.metabolites:
    if '_LLA' in met.id:
        met.id = met.id.replace('_LLA', '_LRE')
        met.name = met.name.replace('_LLA', '_LRE')


'''


print(len(mymodel.reactions),len(mymodel.genes))

myrealist = [i.id for i in mymodel.reactions]
addrea = True
addlist = []


if addrea:

    mymodel.objective =  "ATPM"
    solution = gapfill(mymodel, iNF517, demand_reactions=False);

    for rea in solution[0]:
        rea.notes['from'] = ['iNF517','gap']
        if rea.id not in myrealist:
            mymodel.add_reaction(rea)

            addlist.append(rea.id)
            myrealist.append(rea.id)


    iNF517.objective = "BIOMASS_LLA"
    linear_reaction_coefficients(iNF517)
    f = flux_variability_analysis(iNF517)
    nesf  = f[(f['minimum']>0.0001) | (f['maximum']< -0.0001)]
    myrealist = [i.id for i in mymodel.reactions]

    for i in nesf.index:
        #add fvar
        if i not in myrealist:
            rea = iNF517.reactions.get_by_id(i)
            rea.gene_reaction_rule = ''
            rea.notes['from'] = ['iNF517','gap']
            mymodel.add_reaction(rea)
            addlist.append(rea.id)
            myrealist.append(rea.id)


    #add FBA
    pfba_solution = cobra.flux_analysis.pfba(iNF517)
    need_fluxes = pfba_solution.fluxes[abs(pfba_solution.fluxes) > 1e-10]
    for need_id in need_fluxes.index:
        #if need_id not in myrealist:
        rea = iNF517.reactions.get_by_id(need_id)
        rea.gene_reaction_rule = ''
        rea.notes['from'] = ['iNF517', 'gap']
        mymodel.add_reaction(rea)
        addlist.append(need_id)
        myrealist.append(need_id)


    mymodel.objective =  "BIOMASS_LLA"
    solution = mymodel.optimize()

    print(solution)
    print(iNF517.optimize())
    for i in open('Lreuold.txt'):
        str = i
        break
    addlist2 = []
    for i in addlist:
        if i in str:
            rea = mymodel.reactions.get_by_id(i)
            rea.notes['from'].append('Lreuold')
        else:
            addlist2.append(i)


My_def.io_outtxt(mymodel,'mymodle.txt')




'''

midume_dic = {}
for i in mymodel.reactions:
    if "EX_" in i.id:
        midume_dic[i.id] = [i.bounds]
        
with open('midume.pickle', 'wb') as f:
    pickle.dump(midume_dic, f)



'''











