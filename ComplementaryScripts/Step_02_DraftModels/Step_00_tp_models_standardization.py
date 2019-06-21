#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 2019-05-17

"""Step_00_tp_models_standardization.py
:description : script to process the template models (iBT517 & iNF517 & iML1515)
:param :  template models
:returns:  standardlized models
:rtype: model

fixed cases:
1. iBT721 mets id have 'LSQBKT' or '_RSQBKT' removed
2. all models compared with bigg database if id not in current bigg(Apr 2019) but in old bigg, replace by current bigg id
3. iNF517 reaid have '_copy1' or '_copy_2' substr, compared and keep remove one
4. collect manuallist according reports

"""


import os

import cobra
import re
from My_def.model_report import *
import My_def
import pandas as pd

pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)


def iNF517_report(iNF517, bigg_rea_df, rea_report_df):
    '''
    for iNF517
    reaction id contain _copy1 and _copy2 problems
    :param model:
    :param bigg_rea_df:
    :param rea_report_df:
    :return:
    '''
    model = iNF517.copy()
    columns = list(rea_report_df.columns)
    temp_df = rea_report_df[rea_report_df['descripation'] == 'not in bigg']
    temp_df = temp_df.copy()

    for index in range(len(temp_df)) :
        row_list = list(temp_df.iloc[index])
        id_in_tp = row_list[columns.index('id_in_tp')]

        if '_copy' in  id_in_tp :
            new_id = id_in_tp.split('_copy')[0]
            row_list[columns.index('notes')] = new_id

            rea = model.reactions.get_by_id(id_in_tp)
            rea.id = new_id
            temp_row_list = each_rea_report(rea, model, columns, bigg_rea_df)
            rea.id = id_in_tp

            row_list = merge_row_list(row_list , temp_row_list,columns)
            row_list[columns.index('notes')] = new_id

            if 'same' in row_list[columns.index('descripation')]:
                row_list[columns.index('notes')] = 'keep'

            temp_df.iloc[index] = row_list
            #print(temp_df.iloc[index])

            #rea_report_df.loc[rea_report_df['id_in_tp'] == id_in_tp] = temp_df.iloc[index]

    rea_report_df.update(temp_df)

    # tp_model = cobra.Model('tp')
    # for rea in model.reactions:
    #     if rea.id.endswith('_1'):
    #         tp_model.add_reaction(rea)
    #         tp_model.reactions.get_by_id(rea.id).id = re.sub('_1$', '', rea.id)
    #
    # rea_report_df_2 = rea_report(tp_model, bigg_rea_df)
    # rea_report_df_2['check'] = 'changed_id'
    # rea_report_df.update(rea_report_df_2)

    return rea_report_df

def iNF517_process(iNF517,rea_report_df):

    '''
    chaneg model accordng to the report:
    combine _copy1 and _copy2
    if rea == rea in bigg: keep
    elif bounds != (0,0) :keep
    elif copy1 == copy2 merge gpr,keep one

    chaneg id without _copy and _1
    '''

    df = rea_report_df[rea_report_df['notes'] != '']
    df = df.copy()

    realist = list(df['id_in_tp'])
    for index in range(0, len(realist), 2):
        keep1 = True
        rea_copy1 = iNF517.reactions.get_by_id(realist[index])
        rea_copy2 = iNF517.reactions.get_by_id(realist[index + 1])
        # if 'same' in df.iloc[index]['descripation']:
        #     rea = iNF517.reactions.get_by_id(realist[index])
        #     rea_copy = iNF517.reactions.get_by_id(realist[index + 1])
        # elif 'same' in df.iloc[index + 1]['descripation']:
        #     keep1 = False
        #     rea = iNF517.reactions.get_by_id(realist[index + 1])
        #     rea_copy = iNF517.reactions.get_by_id(realist[index])


        if rea_copy1.bounds == (0.0, 0.0) or (rea_copy1.bounds[0] >= rea_copy2.bounds[0] and rea_copy1.bounds[1] <= rea_copy2.bounds[1]) :
            keep1 = False
            rea = rea_copy2
            rea_copy = rea_copy1
        elif rea_copy2.bounds == (0.0, 0.0):
            rea = rea_copy1
            rea_copy = rea_copy2
        else:

            rea = rea_copy1
            rea_copy = rea_copy2


        if rea.gene_reaction_rule == '':
            rea.gene_reaction_rule = rea_copy.gene_reaction_rule
        elif rea_copy.gene_reaction_rule == '':
            pass
        elif rea.gene_reaction_rule != rea_copy.gene_reaction_rule:
            rea.gene_reaction_rule = '(' + rea.gene_reaction_rule + ') or (' + rea_copy.gene_reaction_rule + ')'

        if keep1:
            df['notes'].iloc[index] = 'keep'
            df['notes'].iloc[index + 1] = 'replaced'
        else:
            df['notes'].iloc[index + 1] = 'keep'
            df['notes'].iloc[index] = 'replaced'

        rea.id = realist[index].split('_copy')[0]
        df['new_id'].iloc[index + 1] = rea.id
        df['new_id'].iloc[index] = rea.id

        rea_copy.remove_from_model()
        #print(index)

    rea_report_df.update(df)
    return iNF517


if  __name__ == '__main__':
    os.chdir('../../ComplementaryData/Step_02_DraftModels/Template/template_models/')
    os.system('cp -r ../../../Initial_data/template_models/ ./')

    bigg_rea_df = pd.read_csv('../../../bigg_database/bigg_rea_df.csv', sep='\t')
    bigg_met_df = pd.read_csv('../../../bigg_database/bigg_met_df.csv', sep='\t')
    #bigg_model = cobra.io.read_sbml_model('../../../bigg_database/universe_draft.xml')

    case = ['iBT721','iNF517','iML1515']
    # %% <iBT721 standardization>
    if 'iBT721' in case:
        print('----- Processing iBT721 -----')
        iBT721 = cobra.io.read_sbml_model('iBT721.xml')
        iBT721.id = 'iBT721'

        # pre processing case 'LSQBKT' in rea id
        for met in iBT721.metabolites:
            met.id = met.id.replace('LSQBKT', '')
            met.id = met.id.replace('_RSQBKT', '')

        # integrate the ids
        met_report_df = met_report(iBT721, bigg_met_df, compartment='_')
        iBT721 = standardlize_met(iBT721, met_report_df)
        rea_report_df = rea_report(iBT721, bigg_rea_df)
        iBT721 = standardlize_rea(iBT721, rea_report_df)

        # %% save report
        iBT721_initial_report = met_report_df.append(rea_report_df)
        iBT721_initial_report.to_csv('iBT721_initial_report.csv', sep='\t', index=False)
        #cobra.io.save_json_model(iBT721, iBT721.id + '_standlized.json')

    # %% <iNF517 standardization>
    if 'iNF517' in case:
        print('----- Processing iNF517 -----')
        iNF517 = cobra.io.read_sbml_model('iNF517.xml')
        iNF517.id = 'iNF517'

        met_report_df = met_report(iNF517, bigg_met_df, compartment='_')
        iNF517 = standardlize_met(iNF517, met_report_df)

        rea_report_df = rea_report(iNF517, bigg_rea_df)
        rea_report_df = iNF517_report(iNF517, bigg_rea_df, rea_report_df)

        iNF517 = iNF517_process(iNF517, rea_report_df)
        # TODO: cheeck reactions that have same id but different mets
        iNF517 = standardlize_rea(iNF517, rea_report_df)

        # %% <step change iNF517 special mets id >
        #case _LLA sepcial mets(biomass composition) or id
        for rea in iNF517.reactions:
            if '_LLA' in rea.id:
                rea.id = rea.id.replace('_LLA','_LRE')
                rea.name = rea.name.replace('_LLA','_LRE')
        for met in iNF517.metabolites:
            if '_LLA' in met.id:
                met.id = met.id.replace('_LLA', '_LRE')
                met.name = met.name.replace('_LLA', '_LRE')

        # %% save
        iNF517_initial_report = met_report_df.append(rea_report_df)
        iNF517_initial_report.to_csv('iNF517_initial_report.csv', sep='\t', index=False)
        #cobra.io.save_json_model(iNF517, iNF517.id + '_standlized.json')

    # %% <iML1515 standardization>
    if 'iML1515' in case:
        print('----- Processing iML1515 -----')
        iML1515 = cobra.io.read_sbml_model('iML1515.xml')
        iML1515.id = 'iML1515'

        met_report_df = met_report(iML1515, bigg_met_df, compartment='_')
        iML1515 = standardlize_met(iML1515, met_report_df)
        rea_report_df = rea_report(iML1515, bigg_rea_df)
        iML1515 = standardlize_rea(iML1515, rea_report_df)

        # save
        iML1515_initial_report = met_report_df.append(rea_report_df)
        iML1515_initial_report.to_csv('iML1515_initial_report.csv', sep='\t', index=False)
        #cobra.io.save_json_model(iML1515, iML1515.id + '_standlized.json')

    # %% <Manual handling> Manual handling according to the reports
    #
    print('----- Manual handling -----')
    modellist = ['iNF517','iBT721','iML1515']

    # remove met id duplactions
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
            locals()[model] = My_def.merge_model.merge_metabolitesid(locals()[model], k, v)
    # %% save
    print('----- Save models -----')
    cobra.io.save_json_model(iBT721, iBT721.id + '_standlized.json')
    cobra.io.save_json_model(iNF517, iNF517.id + '_standlized.json')
    cobra.io.save_json_model(iML1515, iML1515.id + '_standlized.json')

    # %% <seq processing>

    print('----- seq processing -----')
    os.chdir('../template_seqs/')

    # %%
    os.system('cp -r ../../../Initial_data/template_seqs/ ./')

    templatelist = ['iBT721',
                    'iNF517',
                    'iMP429',
                    'iYO844',
                    'iML1515']

    for i in templatelist:
        gbk_file = i + '.gbff'
        faa_file = i + '.faa'
        My_def.seq_ana.gbk2faa(gbk_file,faa_file,locus_tag = 'locus_tag')


    print('===== Template models and seqs processed, Done =====')






