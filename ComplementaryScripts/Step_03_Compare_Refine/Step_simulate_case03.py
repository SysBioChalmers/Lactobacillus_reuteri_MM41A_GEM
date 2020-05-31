#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 2019-09-02

"""Step_simulate.py
:description : script
:param : 
:returns: 
:rtype: 
"""

import os

import cobra
import matplotlib.pyplot as plt
import numpy as np

os.chdir('../../ComplementaryData/Step_03_Compare_Refine/')
print('----- loading model -----')
iHL622 = cobra.io.load_json_model('../../ModelFiles/iHL622.json')

# %% <production ability >
print('----- change medium -----')
iHL622.objective = "BIOMASS"
experiment_medium = {
    'EX_glc__D_e': [-20, -20, ],
    'EX_glyc_e': [0.00, -20, ]}

predict_result = {'EX_lac__L_e': [],
                  'EX_ac_e': [],
                  'EX_etoh_e': [],
                  'EX_13ppd_e': [],
                  }

for i in range(0, len(experiment_medium['EX_glc__D_e'])):#len(experiment_result)
    model = iHL622.copy()

    for rea in experiment_medium.keys():        # close other carbon sourse
        bound = experiment_medium[rea][i]
        model.reactions.get_by_id('EX_cys__L_e').bounds = (0, 1000)
        model.reactions.get_by_id('EX_gly_e').bounds = (0, 1000)
        model.reactions.get_by_id('EX_ala__L_e').bounds = (0, 1000)
        model.reactions.get_by_id('EX_leu__L_e').bounds = (0, 1000)
        model.reactions.get_by_id('EX_ile__L_e').bounds = (0, 1000)
        model.reactions.get_by_id('EX_thr__L_e').bounds = (0, 1000)
        model.reactions.get_by_id('EX_arg__L_e').bounds = (0, 1000)
        model.reactions.get_by_id('EX_asn__L_e').bounds = (0, 1000)
        model.reactions.get_by_id('EX_asp__L_e').bounds = (0, 1000)
        model.reactions.get_by_id('EX_glu__L_e').bounds = (0, 1000)
        model.reactions.get_by_id('EX_gln__L_e').bounds = (0, 1000)
        model.reactions.get_by_id('EX_glu__L_e').bounds = (0, 1000)
        model.reactions.get_by_id('EX_ser__L_e').bounds = (0, 1000)
        model.reactions.get_by_id('EX_his__L_e').bounds = (0, 1000)

        if bound <= 0:
            model.reactions.get_by_id(rea).bounds = (bound, 0)
        else:
            model.reactions.get_by_id(rea).bounds = (0, bound)
    for boj in predict_result.keys():
        model.objective = boj
        sol = model.optimize()
        prod_i = model.metabolites.get_by_id(boj.replace('EX_', '')).formula_weight
        # if boj != 'EX_13ppd_e':
        #     sub = experiment_medium['EX_glc__D_e'][i] * 108
        # else:
        #     sub = experiment_medium['EX_glyc_e'][i] * 75
        sub = experiment_medium['EX_glc__D_e'][i] * (-108)
        predict_result[boj].append(round(sol.objective_value * prod_i / sub, 3))
        # predict_result[boj].append(round(sol.objective_value, 3))
# print('Experiment Biomass:', experiment_result)
print('iHL622 :', predict_result)

# %%initial experimentdata:
data_1 = {'EX_glc__D_e': [-17.65, -24.09, -22.23,  -19.27, -17.96],
    'EX_glyc_e': [0.00, 0.00, 0.00,  -21.43, -18.94],
    'EX_lac__L_e': [17.58, 24.08, 20.86, 18.16, 16.97],
    'EX_ac_e': [-0.19, 0.23, 0.44,  9.29, 8.43],
    'EX_etoh_e': [14.48, 20.05, 16.64,  6.25, 6.46],
    'EX_13ppd_e': [0.00, 0.00, 0.00,  17.50, 15.62],}

for key in data_1.keys():
    data_1[key] = np.array(data_1[key])*iHL622.metabolites.get_by_id(key.replace('EX_','')).formula_weight
    if key not in ['EX_glc__D_e','EX_glyc_e']:
        temp = data_1[key]/data_1['EX_glc__D_e']
        data_1[key] = -temp
print('Experiment :', data_1)

data_2 = {
    'EX_lac__L_e': [18.215,17.881],
    'EX_ac_e': [17.058, 19.7],
    'EX_etoh_e': [5.135, 2.558],
    'EX_13ppd_e': [0.00, 0.00 ]}
data_2_err = {
    'EX_lac__L_e': [0.0115,0.095],
    'EX_ac_e': [0.039,0.067],
    'EX_etoh_e': [0.022,0.011],
    'EX_13ppd_e': [0.00, 0.00],}
for key in data_2.keys():
    data_2[key] = np.array(data_2[key])/20
    data_2_err[key] = np.array(data_2_err[key])/20
print('Experiment :', data_2)

data_3 = {
    'EX_lac__L_e': [0.709,0.419],
    'EX_ac_e': [0.333, 0.172],
    'EX_etoh_e': [0.088, 0.059],
    'EX_13ppd_e': [0.328, 0.414 ]}

data_3_err = {
    'EX_lac__L_e': [0.002,0.021],
    'EX_ac_e': [0.005,0.001],
    'EX_etoh_e': [0.001,0.003],
    'EX_13ppd_e': [0.028, 0.031],}
print('Experiment :', data_3)

# %% <draw>
import statistics
predict_result_draw = []
predict_result_draw_gly = []
for k,v in predict_result.items():
    predict_result_draw.append(v[0])
    predict_result_draw_gly.append(v[1])


import brewer2mpl

fig, ax = plt.subplots(figsize=(6, 4))
ax2 = ax.twinx()
bmap = brewer2mpl.get_map('Set2', 'qualitative', 7)
colors = bmap.mpl_colors

# plt.ylim((0.0, 1.0))
x = np.arange(0, 4)
width = 0.3  # the width of the bars

rects1 = ax.bar(x - width / 2, predict_result_draw, width, label='Glucose',
                 color=colors[1])  #
rects2 = ax.bar(x + width / 2, predict_result_draw_gly, width, label='Glucose+Glycerol', color=colors[0])  # ,
rects3 = ax.errorbar(0,-1,yerr = 0,label='Experiment data',fmt='.k')
keys = ['EX_lac__L_e','EX_ac_e','EX_etoh_e','EX_13ppd_e']
for i in range(0,len(keys)):
    key = keys[i]
    rects4 = ax.errorbar([x[i]-0.1 - width/2], np.mean(data_1[key][0:3]),yerr = statistics.stdev(data_1[key][0:3]), fmt='.k')  # ,
    rects4 = ax.errorbar([x[i]-0.1 + width/2], np.mean(data_1[key][3:5]),yerr = statistics.stdev(data_1[key][3:5]), fmt='.k')  # ,

    rects4 = ax.errorbar([x[i] - width/2], data_2[key][0], yerr = data_2_err[key][0],  fmt='.k')  # ,
    rects4 = ax.errorbar([x[i] + width/2], data_2[key][1], yerr = data_2_err[key][1],  color=colors[2],fmt='.k')  # ,

    # rects3 = ax.errorbar([x[i]+0.1 - width/2]*2, data_3[key],yerr = data_3_err[key],color=colors[2] ,fmt='.k')  # ,
    rects4 = ax.errorbar([x[i]+0.1 + width/2]*2, data_3[key],yerr = data_3_err[key],color=colors[2],fmt='.k')  # ,


# rects1_ = ax2.bar(0, 0, label='Model Growth rate', color=colors[0], )
# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Yield (g/g)', fontsize=16)  # color = 'tab:blue'
# ax.tick_params(axis='y')  # , labelcolor='tab:blue'
ax.set_ylim((0, 2.75))

ax.set_title('Model yield ability of different productions', fontsize=18)
labels =['Lactate','Acetate','Ethanol','1,3-propanediol']
ax.set_xticks(x)
ax.set_xticklabels(labels, fontsize=12)
ax.legend(loc='best', fontsize=11)

fig.tight_layout()
plt.show()
fig.savefig('Growth rate simulation case3.png')
