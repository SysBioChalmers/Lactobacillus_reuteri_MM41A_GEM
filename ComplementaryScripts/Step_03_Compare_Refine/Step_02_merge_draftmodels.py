#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 2019-05-28

"""Step_02_merge_draftmodels.py
:description : script
:param : 
:returns: 
:rtype: 
"""
import cobra
import os
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
    notes = ['']*3
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
    else:
        notes[0] = 'mets different'

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

def merge_draftmodels(model1, model2):
    model = model1.copy()
    reaset = set([i.id for i in model.reactions])
    for rea2 in model2.reactions:
        if rea2.id not in reaset:
            model.add_reaction(rea2)
            reaset.add(rea2.id)
        else:
            rea1 = model1.reactions.get_by_id(rea2.id)
            gpr,notes = compare_euqations(rea1, rea2)
            if notes[0] == 'mets different':
                print(rea1,'\n',rea2,'\n',notes[0],'ignore')
            elif notes[2] == '':
                model.reactions.get_by_id(rea2.id).gene_reaction_rule = gpr
    return model


if __name__=='__main__':
    os.chdir('../../ComplementaryData/Step_02_DraftModels/')


    Lreu_ca_standardlized = cobra.io.load_json_model('CarveMe/Lreu_ca_standardlized.json')
    Lreu_ca_gp_standardlized = cobra.io.load_json_model('CarveMe/Lreu_ca_gp_standardlized.json')
    Lreu_from_iNF517 = cobra.io.load_json_model('Template/Lreu_from_iNF517.json')
    Lreu_from_iBT721 = cobra.io.load_json_model('Template/Lreu_from_iBT721.json')

    # %%    <option 1>
    # by cobra function
    #model1 = Lreu_from_iBT721.merge(Lreu_from_iNF517,inplace=False)
    #model1 = model1.merge(Lreu_ca_gp_standardlized,inplace=False)

    # %%    <option 2>
    # by def function

    model2 = merge_draftmodels(Lreu_from_iBT721,Lreu_from_iNF517)

    model2 = merge_draftmodels(model2,Lreu_ca_gp_standardlized)






    #cobra.io.save_json_model(model,'../../ModelFiles/Lreu.json')
