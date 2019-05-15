#!/usr/bin/env python
#-*- coding: utf-8 -*-
# Created by lhao at 2019-04-21

import cobra
from cobra import Model, Reaction, Metabolite

def merge_model(model1, model2,newmodelname = ''):
    for i in model1.reactions:
        i.notes['from'] = [model1.id]

    for i in model2.reactions:
        i.notes['from'] = [model2.id]

    realist1 = [i.id for i in model1.reactions]
    realist2 = [i.id for i in model2.reactions]

    addlist = set(realist2)-set(realist1)

    newmodel = Model(newmodelname)


    newmodel.add_reactions(model1.reactions)
    for i in addlist:

        newmodel.add_reactions([model2.reactions.get_by_id(i)])

    return newmodel

