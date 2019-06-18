#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 2019-05-28

"""Step_02_merge_draftmodels.py
:description : script to merge Lre_ from iNF517 and iBT721 and iML1515
:param : draft models
:returns: merged model
:rtype: 
"""
import cobra
import os
import pandas as pd
import re
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
from My_def.merge_model import *


if __name__=='__main__':
    os.chdir('../../ComplementaryData/Step_02_DraftModels/')

    # %% load models
    Lreu_ca = cobra.io.load_json_model('CarveMe/Lreu_ca.json')
    Lreu_ca_gp = cobra.io.load_json_model('CarveMe/Lreu_ca_gp.json')
    Lreu_from_iML1515 = cobra.io.load_json_model('Template/Lreu_from_iML1515.json')
    Lreu_from_iNF517 = cobra.io.load_json_model('Template/Lreu_from_iNF517.json')
    Lreu_from_iBT721 = cobra.io.load_json_model('Template/Lreu_from_iBT721.json')
    iNF517 = cobra.io.load_json_model('Template/template_models//iNF517_standlized.json')
    iBT721 = cobra.io.load_json_model('Template/template_models//iBT721_standlized.json')
    iML1515 = cobra.io.load_json_model('Template/template_models//iML1515_standlized.json')

    Lreu_ca.id = 'Lreu_ca'
    Lreu_ca_gp.id = 'Lreu_ca_gp'
    Lreu_from_iML1515.id = 'Lreu_from_iML1515'
    Lreu_from_iNF517.id = 'Lreu_from_iNF517'
    Lreu_from_iBT721.id = 'Lreu_from_iBT721'
    iML1515.id = 'iML1515'
    iNF517.id = 'iNF517'
    iBT721.id = 'iBT721'

    note_model_from(Lreu_ca,Lreu_ca.id)
    note_model_from(Lreu_ca_gp,Lreu_ca_gp.id)
    note_model_from(Lreu_from_iNF517, Lreu_from_iNF517.id)
    note_model_from(Lreu_from_iBT721, Lreu_from_iBT721.id)
    note_model_from(Lreu_from_iML1515, Lreu_from_iML1515.id)

    # %%    <option 1>
    # by cobra function
    #model_1 = Lreu_from_iBT721.merge(Lreu_from_iNF517,inplace=False)
    #model_1 = model1.merge(Lreu_ca_gp_standardlized,inplace=False)

    # %%    <option 2>
    # by def function
    # Note：based on Lreu_from_iNF517 because metabolites notes

    model_2, report_df = merge_draftmodels(Lreu_from_iNF517,Lreu_from_iML1515)
    model_2, report_df1 = merge_draftmodels(model_2, Lreu_from_iBT721)
    #model_2, report_df2  = merge_draftmodels(model_2,Lreu_ca_gp)


    # %% Manual handling change id according to printed reports (report_df and report_df1)

    modellist = ['iNF517','iBT721','iML1515','Lreu_from_iNF517','Lreu_from_iBT721','Lreu_ca','Lreu_ca_gp','model_2','Lreu_from_iML1515']
    #according to report 1
    # replace_list = ['isobut_e','isobut_c', 'isoval_c', '2mpal_c', '3mbal_c', 'orn__L_c', 'btd__RR_c', 'cysth__L_c',
    #                 'ribflvRD_c', '5fothf_c', 'g6p__B_c', 'glcn__D_c', 'MCOOH_c', 'orn__L_c', 'g6p__B_c',
    #                 'dtdp6dm_c', 'dtdpglc_c', 'ugmd_c', 'btd__RR_c', 'isoval_e', 'isobut_e', ]
    # keep_list = ['acgal_e', '2ahhmd_c', 'acgal_e']

    manual_dic = {'g6p_B_c':'g6p__B_c',
                  'g1p_B_c':'g1p__B_c',
                  '5fthf_c':'5fothf_c',
                  'glcn_c':'glcn__D_c',
                  'orn_c':'orn__L_c',
                  'dtdprmn_c':'dtdp6dm_c',
                  'dtdpglu_c':'dtdpglc_c',
                  'mpt_c':'MPT_c',
                  'acgal_e':'acgala_e',
                  'cyst__L_c':'cysth__L_c',
                  'moco_c':'Moco_c',
                  'btd_RR_c':'btd__RR_c',
                  '2mpa_c':'isobut_c',
                  '2mpa_e':'isobut_e',
                  '3mba_e':'isoval_e',
                  'rbflvrd_c':'ribflvRD_c'
                  }

    for model in modellist:
        for k,v in manual_dic.items():
            locals()[model] = merge_metabolitesid(locals()[model], k, v)


    # %% check again

    model_2, report_df = merge_draftmodels(Lreu_from_iNF517,Lreu_from_iML1515)
    model_2, report_df1 = merge_draftmodels(model_2, Lreu_from_iBT721)


    #case '_1' or '_2' in reaid compare and keep only one

    for rea in model_2.reactions:

        if rea.id.endswith('_1') or rea.id.endswith('_2'):
            try:
                rea2 = model_2.reactions.get_by_id(re.sub('_.$','',rea.id))
                try:
                    rea3 = iNF517.reactions.get_by_id(re.sub('_.$','',rea.id))
                    rea4  =iML1515.reactions.get_by_id(re.sub('_.$','',rea.id))
                    print(rea3)
                    print(rea4)
                except :
                    pass
                print(rea)
                print(rea2)


                rea2.reaction = rea.reaction
                rea2.gene_reaction_rule = merge_gprule(rea.gene_reaction_rule, rea2.gene_reaction_rule)
                rea2.notes['from'] = rea2.notes['from']+rea.notes['from']
                rea2.notes['from'] = list(set(rea2.notes['from']))
                if rea.bounds!=(0,0):
                    rea2.bounds = rea.bounds
                rea.remove_from_model()

                iNF517.reactions.get_by_id(rea.id).remove_from_model()
                try:
                    iNF517.reactions.get_by_id(rea2.id).remove_from_model()
                except :
                    pass
                iNF517.add_reaction(rea2)
            except  KeyError:
                pass

    #case _LLA sepcial mets(biomass composition) or id
    for model in [iNF517,model_2,Lreu_from_iNF517]:
        for rea in model.reactions:
            if '_LLA' in rea.id:
                rea.id = rea.id.replace('_LLA','_LRE')
                rea.name = rea.name.replace('_LLA','_LRE')
        for met in model.metabolites:
            if '_LLA' in met.id:
                met.id = met.id.replace('_LLA', '_LRE')
                met.name = met.name.replace('_LLA', '_LRE')

# %% change someting wrong:
    # CBMKr_copy1(0, 1000) and CBMKr_copy2(-1000,1000)  bounds different and removed wrong one
    iNF517.reactions.get_by_id('CBMKr').lower_bound = -1000
    Lreu_from_iNF517.reactions.get_by_id('CBMKr').lower_bound = -1000
    model_2.reactions.get_by_id('CBMKr').lower_bound = -1000

    # standlized iNF517 is the same as initial iNF517
    iNF517_initial = cobra.io.read_sbml_model('/Users/lhao/Box Sync/Projects/Project_Lreuteri/Lactobacillus_reuteri_MM41A_GEM/ComplementaryData/Step_02_DraftModels/Template/template_models/iNF517.xml')
    iNF517_initial.objective = 'BIOMASS_LLA'
    iNF517.objective = 'BIOMASS_LRE'

    solution1 = iNF517.optimize()
    solution2 = iNF517_initial.optimize()
    if solution1.objective_value == solution2.objective_value:
        print('same')
    else:
        print('check')

# %% check duplaction

    rea_id = []
    rea_mets = []

    for rea in model_2.reactions:
        rea_id.append(rea.id)
        rea_me  =set()
        for  i in rea.metabolites:
            rea_me.add(i.id)
        rea_mets.append(str(rea_me))

    rea_pd = pd.DataFrame(zip(rea_id,rea_mets),columns=['id','mets']).sort_values(by = ['mets'])

    check_id = rea_pd[rea_pd['mets'].duplicated(keep = False)]['id']
    check_id = list(check_id)

    #Manual Check:
    for i in list(check_id):
        rea = model_2.reactions.get_by_id(i)
        print(rea,rea.bounds,rea.gene_reaction_rule)
    keep_move_dic = {'APATi':'APAT',
                 'GNKr':'GNK',
                 'FOLD3':'DHPS3',
                 'MPL':'MLTP4',
                 'TDPGDH':'TDPGDH_1',
                 'G3PD1ir':'G3PD1',
                 'P5CR':'P5CRr',
                 'PROTRS_1':'PROTRS',
                 'TYRTRS_1':'TYRTRS',
                 'PSUDS':'YUMPS'}

    for k,v in keep_move_dic.items():
        rea1 = model_2.reactions.get_by_id(k)
        rea2 = model_2.reactions.get_by_id(v)
        rea1.gene_reaction_rule = merge_gprule(rea1.gene_reaction_rule, rea2.gene_reaction_rule)
        rea1.notes['from'] = list(set(rea1.notes['from']+rea2.notes['from']))
        model_2.reactions.get_by_id(v).remove_from_model()

    # case _p compartment removed
    # remove —— reactions
    pset = set()

    for i in model_2.metabolites:
        if i.id.endswith('_p'):
            for rea in i.reactions:
                print(rea)
                pset.add(rea.id)
    for i in pset:
        model_2.reactions.get_by_id(i).remove_from_model()


    for met in model_2.metabolites:
        if len(met.reactions) ==0:
            print('onemet',met)
            met.remove_from_model()

    removegenlist = []
    for gene in model_2.genes:
        if len(gene.reactions) ==0:
            print('onegene',gene)
            removegenlist.append(gene)

    cobra.manipulation.remove_genes(model_2,removegenlist)


    # for i in range(0,len(a),2):
    #     rea1 = model_2.reactions.get_by_id(a[i])
    #     rea2 = model_2.reactions.get_by_id(a[i+1])
    #     if rea1.metabolites == rea2.metabolites:
    #
    #         if rea1.bounds == rea2.bounds:
    #             bounds = rea2.bounds
    #         else:
    #             lbound = min(rea1.bounds[0],rea2.bounds[0])
    #             ubound = max(rea1.bounds[1],rea2.bounds[1])
    #             bounds = (lbound,ubound)
    #
    #         gpr = merge_gprule(rea1.gene_reaction_rule, rea2.gene_reaction_rule)
    #         from_note = list(set(rea1.notes['from'].expend(rea2.notes['from'])))
    #
    #         if '_' in rea1.id:
    #             remove = 1
    #         else:
    #             remove = 2
    #
    #
    #     else:
    #         print('check', rea1)
    #         print('check', rea2)
    #
    #
    #
    #     rea = model_2.reactions.get_by_id(a[i])
    #     print(rea,rea.bounds,rea.gene_reaction_rule)








    #%% save models
    cobra.io.save_json_model(model_2,'../Step_03_Compare_Refine/Lreu_merged.json')
    cobra.io.save_json_model(Lreu_from_iML1515,'../Step_03_Compare_Refine/Lreu_from_iML1515.json')
    cobra.io.save_json_model(Lreu_from_iNF517,'../Step_03_Compare_Refine/Lreu_from_iNF517.json')
    cobra.io.save_json_model(Lreu_from_iBT721,'../Step_03_Compare_Refine/Lreu_from_iBT721.json')
    cobra.io.save_json_model(iNF517,'../Step_03_Compare_Refine/iNF517.json')
    cobra.io.save_json_model(iBT721,'../Step_03_Compare_Refine/iBT721.json')
    cobra.io.save_json_model(iML1515,'../Step_03_Compare_Refine/iML1515.json')



