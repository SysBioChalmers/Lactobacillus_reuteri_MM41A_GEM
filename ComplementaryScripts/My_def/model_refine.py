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
    for met in model.metabolites:
        if len(met.reactions) ==0:
            #print('onemet',met)
            met.remove_from_model()
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

