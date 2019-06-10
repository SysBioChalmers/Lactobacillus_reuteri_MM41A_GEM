#!/usr/bin/env python
#-*- coding: utf-8 -*-
# Created by lhao at 2019-04-21

import cobra
from cobra import Model, Reaction, Metabolite
import os
import pandas as pd
import re
import imp

def dir_rea_metabolites(rea_metbolites):
    '''
    convert reaction.metabolites to str

    :param rea_metbolites: reaction.metabolites a special dictionary
    :return: a str of rea_metbolites
    eval(str_rea_metabolites(rea_metbolites)) could be a normal dictionary
    '''
    temp_dic = {}
    for k, v in rea_metbolites.items():
        temp_dic[k.id] = v
    return temp_dic

def dir_rea_revers_metabolites(rea_metbolites):

    temp_dic = {}
    for k, v in rea_metbolites.items():
        temp_dic[k.id] = -v
    return temp_dic

def compare_euqations(rea1,rea2):
    gpr = rea1.gene_reaction_rule
    notes = ['mets different','bounds different','gpr different']

    if dir_rea_metabolites(rea1.metabolites) == dir_rea_metabolites(rea2.metabolites) :
        notes[0] = 'mets same'
        if rea1.bounds == rea2.bounds:
            notes[1] = 'bounds same'
        if rea1.gene_reaction_rule == rea2.gene_reaction_rule:
            notes[2] = 'gpr same'
        else:
            gpr = merge_gprule(rea1.gene_reaction_rule, rea2.gene_reaction_rule)
    elif dir_rea_metabolites(rea1.metabolites) == dir_rea_revers_metabolites(rea2.metabolites):
        notes[0] = 'mets same'
        if rea1.upper_bound == -rea2.lower_bound and rea1.lower_bound == -rea2.upper_bound:
            notes[1] = 'bounds same'
        if rea1.gene_reaction_rule == rea2.gene_reaction_rule:
            notes[2] = 'gpr same'
        else:
            gpr = merge_gprule(rea1.gene_reaction_rule, rea2.gene_reaction_rule)

    return gpr,notes

def merge_gprule(gpr1,gpr2):
    genset1 = cobra.core.gene.parse_gpr(gpr1)[1]
    genset2 = cobra.core.gene.parse_gpr(gpr2)[1]
    if genset2.issubset(genset1):
        gpr = gpr1
    elif genset1.issubset(genset2):
        gpr = gpr2
    else:
        gpr = '( %s ) or ( %s )'% (gpr1,gpr2)
    return gpr

def merge_draftmodels(model1, model2, inplace = False):

    if inplace:
        model = model1
    else:
        model = model1.copy()

    #model = model1
    reaset = set([i.id for i in model.reactions])
    data = []

    for rea2 in model2.reactions:
        rowlist = ['']*6
        rowlist[0] = rea2.id
        rowlist[1] = 'add'
        if rea2.id not in reaset:
            model.add_reaction(rea2)
            reaset.add(rea2.id)
            rowlist[2] = 'new'
        else:
            rea1 = model1.reactions.get_by_id(rea2.id)
            gpr,notes = compare_euqations(rea1, rea2)
            rowlist[2] = ','.join(notes)

            if notes[0] == 'mets different':
                rowlist[1] = 'skip'

                rowlist[3] = rea1.reaction
                rowlist[4] = rea2.reaction

                dif = [i for i in dir_rea_metabolites(rea1.metabolites).keys() if i not in dir_rea_metabolites(rea2.metabolites).keys()]
                if len(dif)==0:
                    dif = 'coefficient dif'
                rowlist[5] = dif

            elif notes[2] == 'gpr different':
                rowlist[3] = rea1.gene_reaction_rule
                rowlist[4] = rea2.gene_reaction_rule
                model.reactions.get_by_id(rea2.id).gene_reaction_rule = gpr
                model.reactions.get_by_id(rea2.id).notes['from'].extend(rea2.notes['from'])
                model.reactions.get_by_id(rea2.id).notes['from'] = list(set(model.reactions.get_by_id(rea2.id).notes['from']))

            elif notes[1] == 'bounds different':
                rowlist[3] = str(rea1.bounds)
                rowlist[4] = str(rea2.bounds)
                model.reactions.get_by_id(rea2.id).notes['from'].extend(rea2.notes['from'])
                model.reactions.get_by_id(rea2.id).notes['from'] = list(set(model.reactions.get_by_id(rea2.id).notes['from']))
        data.append(rowlist)
    report_df = pd.DataFrame(data,index=None,columns = ['rea_id', 'add_skip','describ','fea1','fea2','dif'])

    print('\033[0;31;47m')
    print('mets different reaction list')
    print(report_df[report_df['describ'].str.contains('mets different')])
    #display(report_df[report_df['describ'].str.contains('mets different')])

    return model.copy(),report_df

def merge_metabolitesid(model1, new_id, old_id):
    model = model1.copy()
    try:
        model.metabolites.get_by_id(old_id).id = new_id

    except ValueError:

        for rea in model.metabolites.get_by_id(old_id).reactions:
            # rea.reaction = re.sub(r'(^|\b)'+id_in_tp+'(\b|$)', id_in_bigg, rea.reaction)
            rea_equ = rea.reaction
            model.reactions.get_by_id(rea.id).reaction = re.sub(r'\b%s\b'%old_id, new_id, rea.reaction)

        model.metabolites.get_by_id(old_id).remove_from_model()
        print(new_id + ' already in model!!!')

    except KeyError:
        print(old_id + ' not in model!!! skiped')

    return model

def judge(model1,old_id,model2,new_id):

    try:
        print(model1.metabolites.get_by_id(old_id).name)
        print(model2.metabolites.get_by_id(new_id).name)

        if model1.metabolites.get_by_id(old_id).formula == '' or model2.metabolites.get_by_id(new_id).formula=='':
            print('Manual handling!!! no formula')
        elif model1.metabolites.get_by_id(old_id).formula == model2.metabolites.get_by_id(new_id).formula:
            #merge_metabolitesid(model1, new_id, old_id)
            return True
        else:
            print('%s and %s not same'%old_id, new_id)
            return False
    except KeyError:
        print(old_id + 'not in model!!! skiped')
        return False


def note_rea_from(rea,notes):
    if 'from' not in rea.notes:
        rea.notes['from'] = [notes]
    else:
        if type(rea.notes['from']) == set:
            rea.notes['from'] = list(rea.notes['from'])
        elif type(rea.notes['from']) == str:
            rea.notes['from'] = [rea.notes['from']]
        rea.notes['from'].append(notes)

    rea.notes['from'] = list(set(rea.notes['from']))

def note_model_from(model,notes):
    for rea in model.reactions:
        rea = note_rea_from(rea, notes)

if __name__=='__main__':
    os.chdir('../../ComplementaryData/Step_02_DraftModels/')

    # %% load models
    Lreu_ca = cobra.io.load_json_model('CarveMe/Lreu_ca.json')
    Lreu_ca_gp = cobra.io.load_json_model('CarveMe/Lreu_ca_gp.json')
    Lreu_from_iNF517 = cobra.io.load_json_model('Template/Lreu_from_iNF517.json')
    Lreu_from_iBT721 = cobra.io.load_json_model('Template/Lreu_from_iBT721.json')

    Lreu_ca.id = 'Lreu_ca'
    Lreu_ca_gp.id = 'Lreu_ca_gp'
    Lreu_from_iNF517.id = 'Lreu_from_iNF517'
    Lreu_from_iBT721.id = 'Lreu_from_iBT721'

    note_model_from(Lreu_ca,Lreu_ca.id)
    note_model_from(Lreu_ca_gp,Lreu_ca_gp.id)
    note_model_from(Lreu_from_iNF517, Lreu_from_iNF517.id)
    note_model_from(Lreu_from_iBT721, Lreu_from_iBT721.id)

    # %%    <option 1>
    # by cobra function
    #model_1 = Lreu_from_iBT721.merge(Lreu_from_iNF517,inplace=False)
    #model_1 = model1.merge(Lreu_ca_gp_standardlized,inplace=False)

    # %%    <option 2>
    # by def function
    # Note：based on Lreu_from_iNF517 because metabolites notes

    #model2,report_df = merge_draftmodels(Lreu_from_iBT721,Lreu_from_iNF517)
    model_2, report_df1 = merge_draftmodels(Lreu_from_iNF517, Lreu_from_iBT721)
    model_2, report_df2  = merge_draftmodels(model_2,Lreu_ca_gp)


    # %% Manual handling

    modellist = ['Lreu_from_iNF517','Lreu_from_iBT721','Lreu_ca','Lreu_ca_gp','model_2']
    #according to report 1

    if judge(model_2, 'g6p_B_c', Lreu_from_iBT721, 'g6p__B_c'):
        for model in modellist:
            locals()[model] = merge_metabolitesid(locals()[model], 'g6p_B_c', 'g6p__B_c')

    if judge(model_2, 'g1p_B_c', Lreu_from_iBT721, 'g1p__B_c'):
        for model in modellist:
            locals()[model] = merge_metabolitesid(locals()[model], 'g1p_B_c', 'g1p__B_c')

    #according to report2
    if judge(model_2, '5fothf_c', Lreu_ca_gp, '5fthf_c'):
        for model in modellist:
            locals()[model] = merge_metabolitesid(locals()[model], '5fthf_c', '5fothf_c')

    if judge(model_2, 'glcn__D_c', Lreu_ca_gp, 'glcn_c'):
        for model in modellist:
            locals()[model] = merge_metabolitesid(locals()[model], 'glcn_c', 'glcn__D_c')

    if judge(model_2, 'orn__L_c', Lreu_ca_gp, 'orn_c'):
        for model in modellist:
            locals()[model] = merge_metabolitesid(locals()[model], 'orn_c', 'orn__L_c')

    if judge(model_2, 'dtdp6dm_c', Lreu_ca_gp, 'dtdprmn_c'):
        for model in modellist:
            locals()[model] = merge_metabolitesid(locals()[model], 'dtdprmn_c', 'dtdp6dm_c')

    #judge(model2, 'glcn__D_e', Lreu_ca_gp, 'dtdprmn_c')

    # %% check again

    model_2, report_df1 = merge_draftmodels(Lreu_from_iNF517, Lreu_from_iBT721)
    model_2, report_df2 = merge_draftmodels(model_2, Lreu_ca_gp)
