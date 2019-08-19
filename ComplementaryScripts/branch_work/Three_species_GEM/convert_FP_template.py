#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 2019-08-14

"""convert_FP_template.py
:description : script
:param : 
:returns: 
:rtype: 
"""

import os
import cobra
from cobra.flux_analysis import flux_variability_analysis
import My_def

os.chdir('../../../ComplementaryData/Initial_data/Three_species_templates')
FP_template = cobra.io.read_sbml_model('iFap484.V01.00.xml')
solution = FP_template.optimize()
print(solution)


solution = flux_variability_analysis(FP_template,fraction_of_optimum=0.9)
nofluxlist  = list(solution[(solution.minimum == 0) and (solution.maximum == 0)].index)

biggmodel = cobra.io.read_sbml_model('../../bigg_database/universe_draft.xml')


model = cobra.Model('FP')

rea_ids = [ rea.id for rea in FP_template.reactions]

biggreaids = [ rea.id for rea in biggmodel.reactions]
# %%

notlist = []
addlist = []



for rea_id in rea_ids:
    if rea_id in biggreaids:
        addlist.append(rea_id)

    else:
        notlist.append(rea_id)

for rea_id in set(addlist) :

    rea = biggmodel.reactions.get_by_id(rea_id)
    model.add_reactions([rea])


for rea_id in set(notlist) - set(nofluxlist):

    rea = FP_template.reactions.get_by_id(rea_id)
    model.add_reactions([rea])


My_def.io_file.model2txt(model,'FP.txt')






