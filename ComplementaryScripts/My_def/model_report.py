#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 2019-05-29

"""model_report.py
:description : script to get a report caompare with bigg database
:param : bigg database
:returns:(dataframe) a report
:rtype: dataframe
"""

import os
import re

import cobra
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib_venn import venn3, venn3_circles
from matplotlib_venn import venn2, venn2_circles
import My_def


# @pysnooper.snoop()
def each_met_report(met_id, bigg_met_df, columns):
    row_list = [''] * len(columns)
    id_in_bigg = met_id
    old_bigg_ids = ''

    if met_id in bigg_met_df['bigg_id'].values:
        description = 'id in bigg'
    else:
        met_pattern = re.compile(r'(^|; )' + met_id + '(;|$)')

        temp_df = bigg_met_df[bigg_met_df['old_bigg_ids'].str.contains(met_pattern, regex=True) == True]

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
    row_list[columns.index('id_in_tp')] = met_id
    row_list[columns.index('id_in_bigg')] = id_in_bigg
    row_list[columns.index('descripation')] = description
    row_list[columns.index('old_bigg_ids')] = old_bigg_ids

    return row_list


def met_report(model, bigg_met_df, compartment=''):
    '''
    get a report to discribe the metabolites features
    descripation : id in bigg or old bigg or not

    :param model:
    :return: met_report_df panda dataframe, a report of metabolites information
    '''

    columns = ['type', 'id_in_tp', 'id_in_bigg', 'descripation', 'old_bigg_ids', 'feature_tp', 'feature_bigg', 'notes']
    met_report_df = pd.DataFrame(columns=columns)

    for met in model.metabolites:
        row_list = each_met_report(met.id, bigg_met_df, columns)

        if row_list[columns.index('descripation')] == 'not in bigg' and compartment != '':
            if compartment in met.id:
                met_id = met.id
                met_id = re.split(compartment + r'.?.$', met_id)[0]
                row_list2 = each_met_report(met_id, bigg_met_df, columns)

                if row_list2[columns.index('descripation')] != 'not in bigg':
                    row_list2[columns.index('descripation')] = 'without compartment(' + row_list2[
                        columns.index('descripation')] + ')'
                    row_list[2:] = row_list2[2:]

        met_report_df.loc[len(met_report_df)] = row_list
        met_report_df = met_report_df.sort_values(by='descripation').reset_index(drop=True)

    return met_report_df


def standardlize_met(model1, met_report_df, keywords=''):
    '''
    change the model according to the report (replace old bigg id )
    :param model: process model,replace id
    :param met_report_df:
    :param keywords:
    :return: model and report
    '''
    model = model1.copy()
    if keywords == '':
        keywords = ['id in old_bigg_ids', 'without compartment(id in bigg)', 'without compartment(id in old_bigg_ids)']

    for keyword in keywords:
        temp_df = met_report_df[met_report_df['descripation'] == keyword]
        temp_df = temp_df.copy()
        for index in range(len(temp_df)):
            id_in_tp = temp_df['id_in_tp'].iloc[index]
            id_in_bigg = temp_df['id_in_bigg'].iloc[index]
            if id_in_bigg != '' and id_in_bigg != id_in_tp:

                try:
                    model.metabolites.get_by_id(id_in_tp).id = id_in_bigg

                except ValueError:

                    for rea in model.metabolites.get_by_id(id_in_tp).reactions:
                        #rea.reaction = re.sub(r'(^|\b)'+id_in_tp+'(\b|$)', id_in_bigg, rea.reaction)
                        model.reactions.get_by_id(rea.id).reaction = re.sub('(^| )'+id_in_tp+'( |$)', id_in_bigg, rea.reaction)
                    model.metabolites.get_by_id(id_in_tp).remove_from_model()
                    print(id_in_bigg + 'already in model!!! combined')
                temp_df['notes'].iloc[index] = 'replaced'

        met_report_df.update(temp_df)
    return model


def str_rea_metabolites(rea_metbolites):
    '''
    convert reaction.metabolites to str

    :param rea_metbolites: reaction.metabolites a special dictionary
    :return: a str of rea_metbolites
    eval(str_rea_metabolites(rea_metbolites)) could be a normal dictionary
    '''
    temp_dic = {}
    for k, v in rea_metbolites.items():
        temp_dic[k.id] = v
    str_mets = str(temp_dic)

    return str_mets


def review_equation(rea, temp_df, row_list, columns):
    '''
    check equation dital incloud euqation, bounds
    :param rea: rea in model
    :param temp_df: a part rea_report that id in bigg/ol_bigg
    :param row_list:
    :param columns:
    :return: row_list to discribe reas same or not

    '''
    rea_met_dic = eval(str_rea_metabolites(rea.metabolites))
    str_bounds = str(rea.bounds)

    if len(temp_df) > 1:

        for i in temp_df['mets']:
            if rea_met_dic == eval(i):
                temp_df = temp_df[temp_df['mets'] == i]
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


# @pysnooper.snoop()
def each_rea_report(rea, model, columns, bigg_rea_df):
    '''
    get a report to discribe the reactions features
    descripation : id in bigg or old bigg , euqation same or not

    :param model:
    :return: rea_report_df panda dataframe, a report of rea information
    :param rea:
    :param columns:
    :param bigg_rea_df:
    :param :
    :return:
    '''
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
    row_list[columns.index('descripation')] = description + row_list[columns.index('descripation')]
    return row_list


# @pysnooper.snoop()
def rea_report(model, bigg_rea_df):
    '''
    get a report to discribe the reactions features
    descripation : id in bigg or old bigg , euqation same or not

    :param model:
    :return: rea_report_df panda dataframe, a report of rea information
    :param rea:
    '''

    columns = ['type', 'id_in_tp', 'id_in_bigg', 'descripation', 'old_bigg_ids', 'feature_tp', 'feature_bigg', 'notes']
    rea_report_df = pd.DataFrame(columns=columns)

    for rea in model.reactions:
        row_list = each_rea_report(rea, model, columns, bigg_rea_df)
        rea_report_df.loc[len(rea_report_df)] = row_list

    rea_report_df = rea_report_df.sort_values(by=['descripation', 'id_in_tp']).reset_index(drop=True)

    return rea_report_df


def standardlize_rea(model1, rea_report_df, keywords=''):
    '''

    change the model according to the report (replace old bigg id )
    :param model: process model,replace id
    :param rea_report_df:
    :param keywords:
    :return: model and report
    '''
    model = model1.copy()
    if keywords == '':
        keywords = ['id in old_bigg(same)', 'id in old_bigg(mets different)',
                    'id in old_bigg(bounds different)', '{id_in_bigg repeat}id in old_bigg(same)']

        temp_df = rea_report_df[(rea_report_df['descripation'].isin(keywords) ) & (rea_report_df['id_in_bigg'] != '')]
        #temp_df = rea_report_df[rea_report_df['id_in_bigg'] != '']
        temp_df = temp_df.copy()
        realist = set([i.id for i in model.reactions])
        for idx, row in temp_df.iterrows():
            id_in_bigg = temp_df.loc[idx, 'id_in_bigg']
            id_in_tp = temp_df.loc[idx, 'id_in_tp']
            if id_in_bigg not in realist:
                model.reactions.get_by_id(id_in_tp).id = id_in_bigg


            temp_df.loc[idx, 'note'] = 'replaced'
            realist.add(id_in_bigg)

        rea_report_df.update(temp_df)
    return model


def merge_row_list(row_list, temp_row_list, columns):
    row_list[columns.index('id_in_bigg')] = temp_row_list[columns.index('id_in_bigg')]
    row_list[columns.index('descripation')] = row_list[columns.index('descripation')] + '[' + temp_row_list[
        columns.index('descripation')] + ']'
    row_list[columns.index('old_bigg_ids'):] = temp_row_list[columns.index('old_bigg_ids'):]
    return row_list


def model_report_compare_bigg(model1, bigg_rea_df, bigg_met_df,compartment='', reaplace_id=True):
    model = model1.copy()
    met_report_df = met_report(model, bigg_met_df, compartment)
    if reaplace_id:
        model = standardlize_met(model, met_report_df, keywords='')
    rea_report_df = rea_report(model, bigg_rea_df)
    if reaplace_id:
        model = standardlize_rea(model, rea_report_df, keywords='')

    model_report = met_report_df.append(rea_report_df, ignore_index=True)
    return model,model_report