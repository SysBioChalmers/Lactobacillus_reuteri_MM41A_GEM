#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 2019-06-21

"""model_refine.py
:description : script functions to refine the models
:param : 
:returns: 
:rtype: 
"""
import cobra
import pandas as pd
def remove_compartment(model1,compartment = '_p'):

    model = model1.copy()
    pset = set()
    for i in model.metabolites:
        if i.id.endswith(compartment):
            for rea in i.reactions:
                #print(rea)
                pset.add(rea.id)
    for i in pset:
        try:
            model.reactions.get_by_id(i).remove_from_model()
        except KeyError:
            print('KeyError ',i)

    model = remove_useless_mets(model)
    model = remove_useless_genes(model)

    return model


def remove_useless_mets(model1):
    model = model1.copy()
    remove_met_list = []
    for met in model.metabolites:
        if len(met.reactions) ==0:
            remove_met_list.append(met.id)
            #print('onemet',met)
    for met in remove_met_list:
        model.metabolites.get_by_id(met).remove_from_model()

    return model

def remove_useless_genes(model1):
    model = model1.copy()
    removegenlist = []
    for gene in model.genes:
        if len(gene.reactions) ==0:
            #print('onegene',gene)
            removegenlist.append(gene)
    cobra.manipulation.remove_genes(model,removegenlist)

    return model


def check_duplicate_rea(model1, cofacters_set,remove = False , ):
    model = model1.copy()
    rea_id = []
    rea_mets = []

    for rea in model.reactions:
        rea_id.append(rea.id)
        rea_me = set()
        for  i in rea.metabolites:
            if i.id not in cofacters_set:
                rea_me.add(i.id)
        rea_mets.append(str(rea_me))

    rea_pd = pd.DataFrame(list(zip(rea_id,rea_mets)),columns=['id','mets']).sort_values(by = ['mets'])

    check_df = rea_pd[rea_pd['mets'].duplicated(keep = False)]
    return check_df

def merge_annotation(annotation_1, annotation_2):
    if set(annotation_1.keys()) & set(annotation_2.keys()) == set([]):
        annotation_temp = dict(annotation_1, **annotation_2)
    else:
        annotation_temp = annotation_1
        for key_i in annotation_2.keys():
            if key_i in annotation_temp.keys():
                ann_1 = annotation_temp[key_i]
                ann_2 = annotation_2[key_i]
                if ann_1 == ann_2:
                    pass
                else:
                    if type(ann_1) == str:
                        ann_1 = [ann_1]
                    if type(ann_2) == str:
                        ann_2 = [ann_2]

                    annotation_temp[key_i] = list(set(ann_1 + ann_2))
            else:
                annotation_temp[key_i] = annotation_2[key_i]
    return annotation_temp


def convert_annotation(annotation):
    annotation_temp = {}
    if len(annotation) == 0:
        return annotation_temp
    replace_dic = {
        'KEGG Compound': 'kegg.compound',
        'CHEBI': 'chebi',
        'BioCyc': 'biocyc',
        'MetaNetX (MNX) Chemical': 'metanetx.chemical',
        'SEED Compound': 'seed.compound',
        'Human Metabolome Database': 'hmdb',
        'MetaNetX (MNX) Equation': 'metanetx.reaction',
        'RHEA': 'rhea',
        'KEGG Reaction': 'kegg.reaction',
        'EC Number': 'ec-code',
        'Reactome':'reactome'
    }
    for i in annotation:
        key_i = i[0]
        val_i = i[1].split('/')[-1]
        if key_i in replace_dic.keys():
            key_i = replace_dic[key_i]
        if key_i in annotation_temp.keys():
            annotation_temp[key_i] = annotation_temp[key_i] + [val_i]
        else:
            annotation_temp[key_i] = [val_i]
    return annotation_temp
