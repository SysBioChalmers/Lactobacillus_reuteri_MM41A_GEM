#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 2019-05-17

"""tp_models_standardlization.py
:description : script to process the template models (iBT517 & iNF517)
:param :  template models
:returns:  standardlized models
:rtype: model
"""

import os
import cobra
import pandas as pd
import pickle
import re
import pysnooper

#@pysnooper.snoop()
def met_report(model,bigg_met_df,reaplace_id = True):
    '''

    :param model:
    :return: met_report_df panda dataframe, a report of metabolites information
    '''

    columns = ['type', 'id_in_tp', 'id_in_bigg', 'descripation', 'old_bigg_ids','feature_tp', 'feature_bigg', 'other']
    met_report_df = pd.DataFrame(columns=columns)

    for met in model.metabolites:

        row_list = ['']*len(columns)
        id_in_bigg = met.id
        old_bigg_ids = ''

        if met.id in bigg_met_df['bigg_id'].values:
            description = 'id in bigg'

        else:
            met_pattern = re.compile(r'(^|; )'+met.id+'(;|$)')
            temp_df = bigg_met_df[bigg_met_df['old_bigg_ids'].str.contains(met_pattern, regex=True) ==True]

            if len(temp_df) != 0:
                description = 'id in old_bigg_ids'
                id_in_bigg = temp_df['bigg_id'].iloc[0]
                old_bigg_ids = temp_df['old_bigg_ids'].iloc[0]

                if len(temp_df) > 1:
                    description = 'more id in old_bigg_ids'
            else:
                description = 'not in bigg'
                id_in_bigg = ''

        row_list[columns.index('type')] = 'met'
        row_list[columns.index('id_in_tp')] = met.id
        row_list[columns.index('id_in_bigg')] = id_in_bigg
        row_list[columns.index('descripation')] = description
        row_list[columns.index('old_bigg_ids')] = old_bigg_ids
        met_report_df.loc[len(met_report_df)] = row_list
        met_report_df = met_report_df.sort_values(by='descripation').reset_index(drop=True)

        if reaplace_id:
            if id_in_bigg != '' and id_in_bigg != met.id:
                met.id = id_in_bigg
    return met_report_df

def standardlize_met(model,met_report_df,keywords = ['id in old_bigg_ids']):
    '''

    :param model: process model,replace id
    :param met_report_df:
    :param keywords:
    :return:
    '''

    for keyword in keywords:
        temp_df = met_report_df[met_report_df['descripation'] == keyword]
        for index in range(len(temp_df)):
            id_in_tp = temp_df['id_in_tp'].icol[index]
            id_in_bigg = temp_df['id_in_bigg'].icol[index]
            if id_in_bigg != '' and id_in_bigg != id_in_tp:
                model.metabolites.get_by_id(id_in_tp).id = id_in_bigg





def str_rea_metabolites(rea_metbolites):
    '''

    :param rea_metbolites: reaction.metabolites a dictionary
    :return: a str of rea_metbolites
    '''
    temp_dic = {}
    for k, v in rea_metbolites.items():
        temp_dic[k.id] = v
    str_mets = str(temp_dic)

    return str_mets

def review_equation(rea, temp_df, row_list,columns):
    # check equation
    rea_met_dic = eval(str_rea_metabolites(rea.metabolites))
    str_bounds = str(rea.bounds)

    if len(temp_df) > 1:

        for i in temp_df['mets']:
            if rea_met_dic == eval(i):
                temp_df  = temp_df[temp_df['mets']==i]
        if len(temp_df) > 1:
            print(temp_df)

    feature_tp = ''
    feature_bigg = ''

    if rea_met_dic == eval(temp_df['mets'].iloc[0]) and str_bounds == temp_df['bounds'].iloc[0]:
        description = '(same)'

    elif rea_met_dic == eval(temp_df['mets'].iloc[0]) and str_bounds != temp_df['bounds'].iloc[0]:
        description = '(bounds different)'
        feature_tp = str_bounds
        feature_bigg = temp_df['bounds'].iloc[0]

    else:
        description = '(mets different)'
        feature_tp = rea.reaction
        feature_bigg = temp_df['equaction_string'].iloc[0]

    row_list[columns.index('descripation')] = description
    if feature_tp != '':
        row_list[columns.index('feature_tp')] = feature_tp
    if feature_bigg != '':
        row_list[columns.index('feature_bigg')] = feature_bigg

    return row_list


#@pysnooper.snoop()
def each_rea_report( rea, columns ,bigg_rea_df,reaplace_id):
    row_list = [''] * len(columns)
    id_in_bigg = rea.id

    if rea.id in bigg_rea_df['bigg_id'].values:

        description = 'id in bigg'
        temp_df = bigg_rea_df[bigg_rea_df['bigg_id'] == rea.id]
        row_list = review_equation(rea, temp_df, row_list, columns)

    else:
        rea_pattern = re.compile(r'(^|; )' + rea.id + '(;|$)')
        temp_df = bigg_rea_df[bigg_rea_df['old_bigg_ids'].str.contains(rea_pattern, regex=True)]

        if len(temp_df) >= 1:
            row_list = review_equation(rea, temp_df, row_list, columns)
            description = 'id in old_bigg'
            row_list[columns.index('old_bigg_ids')] = temp_df['old_bigg_ids'].iloc[0]
            id_in_bigg = temp_df['bigg_id'].iloc[0]
        else:
            description = 'not in bigg'
            id_in_bigg = ''

    row_list[columns.index('type')] = 'rea'
    row_list[columns.index('id_in_tp')] = rea.id
    row_list[columns.index('id_in_bigg')] = id_in_bigg


    if id_in_bigg not in ['', rea.id]:
        realist = [i.id for i in model.reactions]
        if id_in_bigg in realist:
            description = '{id_in_bigg repeat}' + description
        else:
            if reaplace_id:
                rea.id = id_in_bigg

    row_list[columns.index('descripation')] = description + row_list[columns.index('descripation')]
    return row_list



#@pysnooper.snoop()
def rea_report(model, bigg_rea_df,reaplace_id = True):
    columns = ['type', 'id_in_tp', 'id_in_bigg', 'descripation', 'old_bigg_ids','feature_tp', 'feature_bigg', 'other']
    rea_report_df = pd.DataFrame(columns=columns)

    for rea in model.reactions:
        row_list = each_rea_report( rea, columns ,bigg_rea_df,reaplace_id)
        rea_report_df.loc[len(rea_report_df)] = row_list

    rea_report_df = rea_report_df.sort_values(by=['descripation','id_in_tp']).reset_index(drop=True)

    return rea_report_df

def standardlize_rea(model,rea_report_df,keywords = ['id in bigg(same)','id in old_bigg(same)','id in old_bigg(mets different)','id in old_bigg(bounds different)','{id_in_bigg repeat}id in old_bigg(same)']):

    for keyword in keywords:
        temp_df = rea_report_df[rea_report_df['descripation'] == keyword]
        for index in range(len(temp_df)):
            id_in_tp = temp_df['id_in_tp'].iloc[index]
            id_in_bigg = temp_df['id_in_bigg'].iloc[index]

            if id_in_bigg not in ['', id_in_tp] :
                realist = [i.id for i in model.reactions]
                if id_in_bigg not in realist:
                    model.reactions.get_by_id(id_in_tp).id = id_in_bigg

def merge_row_list(row_list , temp_row_list,columns):

    row_list[columns.index('id_in_bigg')] = temp_row_list[columns.index('id_in_bigg')]
    row_list[columns.index('descripation')] = row_list[columns.index('descripation')] +'[' + temp_row_list[columns.index('descripation')] + ']'
    row_list[columns.index('old_bigg_ids'):] = temp_row_list[columns.index('old_bigg_ids'):]
    return row_list

def iNF517_report(model,bigg_rea_df,rea_report_df):
    columns = ['type', 'id_in_tp', 'id_in_bigg', 'descripation', 'old_bigg_ids','feature_tp', 'feature_bigg', 'other']
    temp_df = rea_report_df[rea_report_df['descripation'] == 'not in bigg']

    for index in range(len(temp_df)) :
        row_list = list(temp_df.iloc[index])
        id_in_tp = row_list[columns.index('id_in_tp')]

        if '_copy' in  id_in_tp :
            new_id = id_in_tp.split('_copy')[0]
            row_list[columns.index('other')] = new_id

            rea = model.reactions.get_by_id(id_in_tp)
            rea.id = new_id
            temp_row_list = each_rea_report(rea, columns, bigg_rea_df, reaplace_id = False)
            rea.id = id_in_tp

            row_list = merge_row_list(row_list , temp_row_list,columns)
            row_list[columns.index('other')] = new_id

            if 'same' in row_list[columns.index('descripation')]:
                row_list[columns.index('other')] = 'keep'

            temp_df.iloc[index] = row_list
            #print(temp_df.iloc[index])

            #rea_report_df.loc[rea_report_df['id_in_tp'] == id_in_tp] = temp_df.iloc[index]

    rea_report_df.update(temp_df)

    return rea_report_df


def iNF517_process(iNF517,rea_report_df):

    columns = ['type', 'id_in_tp', 'id_in_bigg', 'descripation', 'old_bigg_ids', 'feature_tp', 'feature_bigg', 'other']
    temp_df = rea_report_df[rea_report_df['other'] != '']
    rea_id_list = list(temp_df['id_in_tp'])

    for i in rea_id_list:
        rea = iNF517.reaction
        print()



if  __name__ == '__main__':
    os.chdir('../../ComplementaryData/02_DraftModels/Template/template_models/')

    # with open('../../../bigg_database/universal_model.pickle', 'rb') as f:
    #    bigg_model = pickle.load(f)
    bigg_rea_df = pd.read_csv('../../../bigg_database/bigg_models_reactions_all.txt')
    bigg_met_df = pd.read_csv('../../../bigg_database/bigg_models_metabolites.txt', sep='\t', usecols=[0, 1, 5])
    t_id = ['iBT721', 'iNF517']
    t_id = ['iNF517']

    for model_id in t_id:
        modelfile = model_id+'.xml'
        model  = cobra.io.read_sbml_model(modelfile)
        if model_id == 'iBT721':
            for met in model.metabolites:
                met.id = met.id.replace('LSQBKT', '')
                met.id = met.id.replace('_RSQBKT', '')

        met_report_df = met_report(model, bigg_met_df, reaplace_id=True)
        standardlize_met(model, met_report_df)
        met_report_df.to_csv(model_id +'_met_report.csv')

        rea_report_df = rea_report(model, bigg_rea_df, reaplace_id=True)

        if model_id == 'iNF517' :
            rea_report_df = iNF517_report(model,bigg_rea_df,rea_report_df)

        standardlize_rea(model, rea_report_df)
        rea_report_df.to_csv(model_id +'_rea_report.csv')

        #cobra.io.save_json_model(model,model_id + '_standlized.json')






