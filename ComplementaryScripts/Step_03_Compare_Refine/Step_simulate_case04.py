#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 4/22/19

"""Step_refine_pipeline_part03_amino_acids.py
:description : script to refine amino acids
:param :
:returns:
:rtype:
"""

import os

import cobra
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

os.chdir('../../ComplementaryData/Step_03_Compare_Refine/')
print('----- loading model -----')
iHL622 = cobra.io.load_json_model('../../ModelFiles/iHL622.json')
iNF517 = cobra.io.read_sbml_model('../Initial_data/template_models/iNF517.xml')
LbReueri = cobra.io.read_sbml_model('../Initial_data/template_models/Lreuteri_530.xml')
iBT721 = cobra.io.load_json_model('../Step_02_DraftModels/Template/template_models/iBT721_standlized.json')
iML1515 = cobra.io.read_sbml_model('../Initial_data/template_models/iML1515.xml')

# %% <amino acids>
unessential = [
    'EX_ala__L_e',
    'EX_asp__L_e',
    'EX_cys__L_e',
    'EX_gly_e',
    'EX_ile__L_e',
    'EX_lys__L_e',
    'EX_pro__L_e',
    'EX_ser__L_e', ]

essential = [
    'EX_arg__L_e',
    # 'EX_asn__L_e',
    # 'EX_gln__L_e',
    'EX_glu__L_e',
    'EX_his__L_e',
    'EX_leu__L_e',
    'EX_met__L_e',
    'EX_phe__L_e',
    'EX_thr__L_e',
    'EX_trp__L_e',
    'EX_tyr__L_e',
    'EX_val__L_e']

aarealist = unessential + essential
aadic = {'None': []}
for aa in aarealist:
    aadic[aa] = []
for model_i in [iHL622, iNF517, LbReueri, iBT721, iML1515]:
    model = model_i.copy()
    model.reactions.get_by_id('EX_glc__D_e').bounds = (-5, 1000)
    model.reactions.get_by_id('EX_glyc_e').bounds = (0, 1000)
    for aa in aarealist:
        if model.reactions.get_by_id(aa).lower_bound == 0:
            model.reactions.get_by_id(aa).lower_bound == -1
    # model.objective = 'BIOMASS'
    solution = model.optimize()
    aadic['None'].append(round(solution.objective_value,3))

    for aa in aarealist:
        model_ = model.copy()
        rea = model_.reactions.get_by_id(aa)
        rea.bounds = (0, 1000)
        if aa == 'EX_asp__L_e':
            model_.reactions.get_by_id('EX_asn__L_e').bounds = (0, 1000)  # asn ~ asp
        elif aa == 'EX_glu__L_e':
            model_.reactions.get_by_id('EX_gln__L_e').bounds = (0, 1000)  # gln ~ glu
        # cobra.flux_analysis.pfba(model)
        solution = model_.optimize()
        aadic[aa].append(round(solution.objective_value, 3))
        # solution = cobra.flux_analysis.pfba(model)  # cobra.flux_analysis.pfba(model)
        # My_def.io_file.solution2txt(solution, model, model.id + aa + '_temp_flux.txt')
    # print(aa, aadic[aa])



# %% plt
import brewer2mpl
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns

experiment_data = [1]+[1] * len(unessential) + [0] * len(essential)
iHL622_data = []
iNF517_data = []
LbReueri_data = []
iBT721_data = []
iML1515_data = []
for key in aadic.keys():
    if aadic[key][0] > 0: iHL622_data.append(1)
    if aadic[key][0] == 0: iHL622_data.append(0)

    if aadic[key][1] > 0: iNF517_data.append(1)
    if aadic[key][1] == 0: iNF517_data.append(0)

    if aadic[key][2] > 0: LbReueri_data.append(1)
    if aadic[key][2] == 0: LbReueri_data.append(0)

    if aadic[key][3] > 0: iBT721_data.append(1)
    if aadic[key][3] == 0: iBT721_data.append(0)

    if aadic[key][4] > 0: iML1515_data.append(1)
    if aadic[key][4] == 0: iML1515_data.append(0)

bmap = brewer2mpl.get_map('Set2', 'qualitative', 7)
colors = bmap.mpl_colors
cm = LinearSegmentedColormap.from_list(
    'cp',[colors[1], colors[0]], N=2)
# data = np.array([experiment_data,iHL622_data,iHL622_data]).T

aalist = ['None'] + [i.split('_')[1] for i in aarealist]
aalist[2] = aalist[2] + '&asn'
aalist[10] = aalist[10] + '&gln'
df = pd.DataFrame({'Experiment': experiment_data,
                   'iHL622': iHL622_data,
                   'iNF517': iNF517_data,
                   'LbReueri': LbReueri_data,
                   'iBT721': iBT721_data,
                   'iML1515': iML1515_data,
                   }, index=aalist)

# dft = df.pivot_table(index='experiment_data', columns='iHL622', values= 'LbReuteri',np.median)

sns.set()

fig, ax = plt.subplots(figsize=(4, 6))

im = sns.heatmap(df, linewidths=.1, ax=ax, cmap=cm,
                 cbar=False)

cbar = plt.colorbar(im.collections[0], #orientation='horizontal',
                    fraction=0.046, pad=0.014, shrink=0.5, aspect=10, )

cbar.set_ticks(np.array([0.25, 0.75]))
cbar.set_ticklabels(('no growth','growth' ))


plt.ylabel('Amino acid omitted')
plt.show()
fig.savefig('Growth rate simulation case4a.png')


