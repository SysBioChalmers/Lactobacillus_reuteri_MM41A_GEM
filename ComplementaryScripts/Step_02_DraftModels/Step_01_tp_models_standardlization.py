#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 2019-05-17

"""Step_01_tp_models_standardlization.py
:description : script to process the template models (iBT517 & iNF517)
:param :  template models
:returns:  standardlized models
:rtype: model
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
        if 'same' in df.iloc[index]['descripation']:
            rea = iNF517.reactions.get_by_id(realist[index])
            rea_copy = iNF517.reactions.get_by_id(realist[index + 1])
        elif 'same' in df.iloc[index + 1]['descripation']:
            keep1 = False
            rea = iNF517.reactions.get_by_id(realist[index + 1])
            rea_copy = iNF517.reactions.get_by_id(realist[index])
        else:
            rea_copy1 = iNF517.reactions.get_by_id(realist[index])
            rea_copy2 = iNF517.reactions.get_by_id(realist[index + 1])

            if rea_copy1.bounds == (0.0, 0.0):
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


# %%
def check_list(model, bigg_model, rea_report_df):
    rea_report_df_tp = rea_report_df[rea_report_df['notes'] == 'manual check'].copy()
    columns = ['type', 'id_in_tp', 'id_in_bigg', 'feature_tp', 'feature_bigg','diff', 'name_tp', 'name_bigg', 'des', 'notes']
    manual_df = pd.DataFrame(columns=columns)

    for index, row in rea_report_df_tp.iterrows():
        row_list = [''] * len(manual_df.columns)
        row_list[columns.index('notes')] = 'please check'
        diff = eval(row['diff'])

        # case 1: coefficient
        # Manual check
        if len(diff[0]) == 0 or len(diff[1]) == 0 or len(diff[0]) != len(diff[1]):
            row['notes'] = row['notes'] + '(coefficient different)'
            case = 'case1'

        elif len(diff[0]) == len(diff[1]):
            # check compounds different

            # case 2: compartments
            tplit1 = set([re.sub('_.$', '', i) for i in diff[0]])
            tplit2 = set([re.sub('_.$', '', i) for i in diff[1]])
            if tplit1 == tplit2:
                case = 'case2'
                row['notes'] = row['notes'] + '(compartments different)'
            # case 3: formula

            else:
                case = 'case3'
                try:
                    tplit1 = set([model.metabolites.get_by_id(i).formula for i in diff[0]])
                    tplit2 = set([bigg_model.metabolites.get_by_id(i).formula for i in diff[1]])

                    if '' in tplit1 | tplit2 or 'X' in tplit1 | tplit2:
                        row['notes'] = row['notes'] + '(formula missing)'

                    elif tplit1 == tplit2 and '' not in tplit1:
                        row['notes'] = row['notes'] + '(formula same)'
                    else:
                        row['notes'] = row['notes'] + '(formula different)'

                except KeyError:
                    row['notes'] = row['notes'] + '(formula missing)'

            for index in range(len(diff[0])):
                id_in_tp = list(diff[0])[index]
                id_in_bigg = list(diff[1])[index]
                row_list[columns.index('type')] = 'met'
                row_list[columns.index('id_in_tp')] = id_in_tp
                row_list[columns.index('id_in_bigg')] = id_in_bigg
                try:
                    row_list[columns.index('feature_tp')] = model.metabolites.get_by_id(id_in_tp).formula
                    row_list[columns.index('feature_bigg')] = model.metabolites.get_by_id(id_in_tp).formula
                    row_list[columns.index('name_tp')] = model.metabolites.get_by_id(id_in_tp).name
                    row_list[columns.index('name_bigg')] = bigg_model.metabolites.get_by_id(id_in_tp).name
                except:
                    pass
                row_list[columns.index('des')] = row['notes']

                if '(formula same)' in row['notes']:
                    row_list[columns.index('notes')] = 'Replace'

        if case == 'case1' or case == 'case2':
            row_list[columns.index('type')] = 'rea'
            row_list[columns.index('id_in_tp')] = row['id_in_tp']
            row_list[columns.index('id_in_bigg')] = row['id_in_bigg']
            row_list[columns.index('feature_tp')] = row['feature_tp']
            row_list[columns.index('feature_bigg')] = row['feature_bigg']
            row_list[columns.index('diff')] = row['diff']
            try:
                row_list[columns.index('name_tp')] = model.reactions.get_by_id(row['id_in_tp']).name
                row_list[columns.index('name_bigg')] = bigg_model.reactions.get_by_id(row['id_in_bigg']).name
            except:
                pass
            row_list[columns.index('des')] = row['notes']

        manual_df.loc[len(manual_df)] = row_list

    rea_report_df.update(rea_report_df_tp)
    return manual_df


def manual_process(manual_df, model1):
    model = model1.copy()
    for index, row in manual_df.iterrows():
        if row['notes'] == 'replace':
            if row['type'] == 'met':
                model = My_def.merge_model.merge_metabolitesid(model, row['id_in_bigg'], row['id_in_tp'])
            elif row['type'] == 'rea':
                model.reactions.get_by_id(row['type']).reaction = row['feature_bigg']
        elif row['notes'] == 'keep':
            if row['type'] == 'met':
                My_def.merge_model.merge_metabolitesid(model, row['id_in_tp'], row['id_in_bigg'])
            elif row['type'] == 'rea':
                model.reactions.get_by_id(row['type']).reaction = row['feature_tp']
    return model


# %%
if  __name__ == '__main__':
    os.chdir('../../ComplementaryData/Step_02_DraftModels/Template/template_models/')

    bigg_rea_df = pd.read_csv('../../../bigg_database/bigg_rea_df.csv', sep='\t')
    bigg_met_df = pd.read_csv('../../../bigg_database/bigg_met_df.csv', sep='\t')
    #bigg_model = cobra.io.read_sbml_model('../../../bigg_database/universe_draft.xml')

    case = ['iBT721','iNF517','iML1515']
    # %% iBT721 report
    if 'iBT721' in case:
        iBT721 = cobra.io.read_sbml_model('iBT721.xml')
        iBT721.id = 'iBT721'

        # pre processing
        for met in iBT721.metabolites:
            met.id = met.id.replace('LSQBKT', '')
            met.id = met.id.replace('_RSQBKT', '')

        # integrate the ids
        met_report_df = met_report(iBT721, bigg_met_df, compartment='_')
        iBT721 = standardlize_met(iBT721, met_report_df)
        rea_report_df = rea_report(iBT721, bigg_rea_df)
        iBT721 = standardlize_rea(iBT721, rea_report_df)

        # %% get a report that need manual check according to the bigg model and formula
        # check reactions that have same id but different mets

        # rea_report_df_tp = rea_report_df.copy()
        # manual_df = check_list(iBT721, bigg_model, rea_report_df_tp)

        # %% check again

        # met_report_df2 = met_report(iBT721, bigg_met_df, compartment='_')
        # iBT721 = standardlize_met(iBT721, met_report_df2)
        # rea_report_df2 = rea_report(iBT721, bigg_rea_df)
        # iBT721 = standardlize_rea(iBT721, rea_report_df2)

        # %% save report
        iBT721_initial_report = met_report_df.append(rea_report_df)
        #iBT721_initial_report.to_csv('iBT721_initial_report.csv', sep='\t', index=False)
        #iBT721_stand_report = met_report_df2.combine(rea_report_df2)
        #iBT721_stand_report.to_csv('iBT721_stand_report.csv', sep='\t', index=False)
        #cobra.io.save_json_model(iBT721, iBT721.id + '_standlized.json')

    # %% iNF517 report

    if 'iNF517' in case:


        iNF517 = cobra.io.read_sbml_model('iNF517.xml')
        iNF517.id = 'iNF517'

        met_report_df = met_report(iNF517, bigg_met_df, compartment='_')
        iNF517 = standardlize_met(iNF517, met_report_df)

        rea_report_df = rea_report(iNF517, bigg_rea_df)
        rea_report_df = iNF517_report(iNF517, bigg_rea_df, rea_report_df)
        iNF517 = iNF517_process(iNF517, rea_report_df)
        # TODO: cheeck reactions that have same id but different mets
        iNF517 = standardlize_rea(iNF517, rea_report_df)

        #met_report_df.to_csv(iNF517.id + '_met_report.csv', sep='\t', index=False)
        #rea_report_df.to_csv(iNF517.id + '_rea_report.csv', sep='\t', index=False)
        iNF517_initial_report = met_report_df.append(rea_report_df)
        iNF517_initial_report.to_csv('iNF517_initial_report.csv', sep='\t', index=False)
        cobra.io.save_json_model(iNF517, iNF517.id + '_standlized.json')
    # %% iML1515
    if 'iML1515' in case:
        iML1515 = cobra.io.read_sbml_model('iML1515.xml')
        iML1515.id = 'iML1515'

        met_report_df = met_report(iML1515, bigg_met_df, compartment='_')
        iML1515 = standardlize_met(iML1515, met_report_df)
        rea_report_df = rea_report(iML1515, bigg_rea_df)
        iML1515 = standardlize_rea(iML1515, rea_report_df)

        iML1515_initial_report = met_report_df.append(rea_report_df)
        #iML1515_initial_report.to_csv('iML1515_initial_report.csv', sep='\t', index=False)
        #cobra.io.save_json_model(iML1515, iML1515.id + '_standlized.json')

    # %% Manual handling according to the report

        # replace_list = ['isobut_e','isobut_c', 'isoval_c', '2mpal_c', '3mbal_c', 'orn__L_c', 'btd__RR_c', 'cysth__L_c',
        #                 'ribflvRD_c', '5fothf_c', 'g6p__B_c', 'glcn__D_c', 'MCOOH_c', 'orn__L_c', 'g6p__B_c',
        #                 'dtdp6dm_c', 'dtdpglc_c', 'ugmd_c', 'btd__RR_c', 'isoval_e', 'isobut_e', ]
        # keep_list = ['acgal_e', '2ahhmd_c', 'acgal_e']
        #
        # manual_df = manual_df.copy()
        # manual_df.loc[manual_df['id_in_tp'].isin(replace_list), 'notes'] = 'replace'
        # manual_df.loc[manual_df['id_in_tp'].isin(keep_list), 'notes'] = 'keep'
        # manual_df = manual_df.sort_values(['type', 'id_in_tp'])
        # manual_df = manual_df.drop_duplicates(subset='id_in_tp', keep='first')
        # manual_df.to_csv('manual_df.scv', sep='\t', index=False)

        #iBT721 = manual_process(manual_df, iBT721)


    iBT721_reas = set([i.id for i in iBT721.reactions])
    iBT721_mets = set([i.id for i in iBT721.metabolites])
    iNF517_reas = set([i.id for i in iNF517.reactions])
    iNF517_mets = set([i.id for i in iNF517.metabolites])
    iML1515_reas = set([i.id for i in iML1515.reactions])
    iML1515_mets = set([i.id for i in iML1515.metabolites])
    #model_2_reas = set([i.id for i in model_2.reactions])
    #model_2_mets = set([i.id for i in model_2.metabolites])

    #only_721_mets = (iBT721_mets - (iNF517_mets|iML1515_mets))&model_2_mets
    #only_517_mets = (iNF517_mets - (iBT721_mets|iML1515_mets))&model_2_mets
    #only_721_reas = (iBT721_reas - (iNF517_reas|iML1515_reas))&model_2_reas
    #only_517_reas = (iNF517_reas - (iBT721_reas|iML1515_reas))&model_2_reas





