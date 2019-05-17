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


def str_rea_metbolites(rea_metbolites):
    '''

    :param rea_metbolites: reaction.metabolites a dictionary
    :return: a str of rea_metbolites
    '''
    temp_dic = {}
    for k, v in rea_metbolites.items():
        temp_dic[k.id] = v
    str_mets = str(temp_dic)

    return str_mets

def model_2_dataframe(model):
    df = pd.DataFrame(columns=['type', 'id_in_tp', 'id_in_bigg', 'descripation', 'feature1', 'feature2', 'other'])
    for met in model.metabolites:
        temp_list = ['']*7
        temp_list[0] = 'metabolites'
        temp_list[1] = met.id



    return df

os.chdir('../../ComplementaryData/02_DraftModels/Template/template_models/')




#with open('../../../bigg_database/universal_model.pickle', 'rb') as f:
#    bigg_model = pickle.load(f)
bigg_rea_df = pd.read_csv('../../../bigg_database/bigg_models_reactions_all.txt')
bigg_met_df = pd.read_csv('../../../bigg_database/bigg_models_metabolites.txt', sep='\t', usecols=[0, 1, 5])


iBT721 = cobra.io.read_sbml_model('iBT721.xml')


dif_df = pd.DataFrame(columns=['type', 'id_in_tp', 'id_in_bigg', 'descripation', 'feature1', 'feature2', 'other'])
for met in iBT721.metabolites:
    bigg_met_df['bigg_id'].isin(met.id)
    pass

#if __name__ == '__main__':


    t_id = ['iBT721','iNF517']