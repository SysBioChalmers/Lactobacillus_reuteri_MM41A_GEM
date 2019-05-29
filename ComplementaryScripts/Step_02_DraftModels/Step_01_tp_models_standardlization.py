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

from My_def.model_report import *


def iNF517_report(model,bigg_rea_df,rea_report_df):
    '''
    for iNF517
    reaction id contain _copy1 and _copy2 problems
    :param model:
    :param bigg_rea_df:
    :param rea_report_df:
    :return:
    '''
    columns = ['type', 'id_in_tp', 'id_in_bigg', 'descripation', 'old_bigg_ids','feature_tp', 'feature_bigg', 'notes']
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

    return rea_report_df


def iNF517_process(iNF517,rea_report_df):

    '''
    chaneg model accordng to the report:
    combine _copy1 and _copy2
    if rea == rea in bigg: keep
    elif bounds != (0,0) :keep
    elif copy1 == copy2 merge gpr,keep one

    chaneg id without _copy
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
        rea_copy.remove_from_model()
        #print(index)

    rea_report_df.update(df)
    return iNF517


if  __name__ == '__main__':
    os.chdir('../../ComplementaryData/Step_02_DraftModels/Template/template_models/')

    bigg_rea_df = pd.read_csv('../../../bigg_database/bigg_rea_df.csv', sep='\t')
    bigg_met_df = pd.read_csv('../../../bigg_database/bigg_met_df.csv', sep='\t')

    # %% iBT721 report
    iBT721 = cobra.io.read_sbml_model('iBT721.xml')
    iBT721.id = 'iBT721'
    for met in iBT721.metabolites:
        met.id = met.id.replace('LSQBKT', '')
        met.id = met.id.replace('_RSQBKT', '')

    met_report_df = met_report(iBT721, bigg_met_df, compartment='_')
    iBT721 = standardlize_met(iBT721, met_report_df)

    # Manual change according the report
    iBT721.metabolites.get_by_id('cysth__L_c').id = 'cyst__L_c'

    rea_report_df = rea_report(iBT721, bigg_rea_df)
    # TODO: cheeck reactions that have same id but different mets
    iBT721 = standardlize_rea(iBT721, rea_report_df)

    met_report_df.to_csv(iBT721.id + '_met_report.csv', sep='\t', index=False)
    rea_report_df.to_csv(iBT721.id + '_rea_report.csv', sep='\t', index=False)
    cobra.io.save_json_model(iBT721, iBT721.id + '_standlized.json')

    # %% iNF517 report
    iNF517 = cobra.io.read_sbml_model('iNF517.xml')
    iNF517.id = 'iNF517'

    met_report_df = met_report(iNF517, bigg_met_df, compartment='_')
    iNF517 = standardlize_met(iNF517, met_report_df)

    rea_report_df = rea_report(iNF517, bigg_rea_df)
    rea_report_df = iNF517_report(iNF517, bigg_rea_df, rea_report_df)
    iNF517 = iNF517_process(iNF517, rea_report_df)
    # TODO: cheeck reactions that have same id but different mets
    iNF517 = standardlize_rea(iNF517, rea_report_df)

    met_report_df.to_csv(iNF517.id + '_met_report.csv', sep='\t', index=False)
    rea_report_df.to_csv(iNF517.id + '_rea_report.csv', sep='\t', index=False)
    cobra.io.save_json_model(iNF517, iNF517.id + '_standlized.json')
