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

# %% <biomass vs od >
print('----- change medium -----')
iHL622.objective = "BIOMASS"

experiment_group = ['A', 'B', 'C', 'D', 'E']
experiment_result = [1.38, 1.88, 1.92, 1.92, 1.90]
experiment_result_err = [0.66, 0.35, 0.69, 0.37, 0.47]
experiment_medium = {
    'EX_glc__D_e': [-20, -20, -20, -20, -20, ],
    'EX_glyc_e': [0.00, -5.00, -5.00, -10.00, -10.],
    'EX_fru_e': [0.00, -1.00, -5.00, -1.00, -5.],
    'EX_lac__L_e': [18.215, 21.334, 20.882, 17.881, 16.577],
    'EX_ac_e': [17.058, 18.301, 18.285, 19.703, 19.643],
    'EX_etoh_e': [5.135, 4.623, 4.312, 2.558, 2.230],
}

# for k in experiment_medium.keys():      # g/L --> mM
#     temp = np.array(experiment_medium[k])*1000/iHL622.metabolites.get_by_id(k.replace('EX_','')).formula_weight
#     experiment_medium[k] = temp

predict_result = []

for i in range(0, len(experiment_result)):
    model = iHL622.copy()
    for rea in experiment_medium.keys():
        bound = experiment_medium[rea][i]
        if bound <= 0:
            model.reactions.get_by_id(rea).bounds = (bound, 0)
        elif bound >= 0:
            model.reactions.get_by_id(rea).bounds = (0, bound)
    sol = model.optimize()
    predict_result.append(round(sol.objective_value, 3))

print('Experiment Biomass:', experiment_result)
print('iHL622 Biomass:', predict_result)

# %% <vitmin B12 > NOTE: error

# experiment_medium = {
#     'BIOMASS': predict_result,
#     'EX_glc__D_e': [-20, -20, -20, -20, -20, ],
#     'EX_glyc_e': [0.00, -5.00, -5.00, -10.00, -10.],
#     'EX_fru_e': [0.00, -1.00, -5.00, -1.00, -5.], }
#
# predict_result_b12 = []
# for i in range(0, len(experiment_result)):
#     model = iHL622.copy()
#     rea = cobra.Reaction('EX_adeadocbl_c')
#     model.add_reaction(rea)
#     model.reactions.get_by_id('EX_adeadocbl_c').reaction = 'adeadocbl_c --> '
#     model.objective = 'EX_adeadocbl_c'
#     # model.reactions.get_by_id('EX_ade_e').bounds = (0,0)
#     for rea in experiment_medium.keys():
#         bound = experiment_medium[rea][i]
#         if rea == 'BIOMASS':
#             model.reactions.get_by_id(rea).bounds = (bound, bound)
#
#         elif bound <= 0:
#             model.reactions.get_by_id(rea).bounds = (bound, 0)
#         elif bound >= 0:
#             model.reactions.get_by_id(rea).bounds = (0, bound)
#     predict_result_b12.append(
#         round(model.optimize().objective_value * 1355.365, 3))  # Cobalamin: Molar mass: 1,355.365 g/mol
# print('iHL622 b12:', predict_result_b12)

# %% <draw>

import brewer2mpl

fig, ax = plt.subplots(figsize=(6, 4))
ax2 = ax.twinx()
bmap = brewer2mpl.get_map('Set2', 'qualitative', 7)
colors = bmap.mpl_colors

# plt.ylim((0.0, 1.0))
x = np.arange(0, 5)
width = 0.25  # the width of the bars

rects2 = ax.bar(x + width / 2, predict_result, width, label='Model Growth rate', color=colors[0])  # ,
rects1 = ax2.bar(x - width / 2, experiment_result, width, yerr=experiment_result_err, label='Experiment OD600',
                 color=colors[1])  #
rects1_ = ax2.bar(0, 0, label='Model Growth rate', color=colors[0], )
# Add some text for labels, title and custom x-axis tick labels, etc.
ax2.set_ylabel("OD600", fontsize=16)
ax.set_ylabel('Growth rate (mmol/gDW/h)', fontsize=16)  # color = 'tab:blue'
# ax.tick_params(axis='y')  # , labelcolor='tab:blue'
ax2.set_ylim((0, 3.2))
ax.set_ylim((0, 2.2))

ax.set_title('Growth rate simulation', fontsize=18)
labels = [''] + experiment_group
ax2.set_xticklabels(labels, fontsize=16)

ax2.legend(loc='best', fontsize=11)
# ax2.legend(loc='best', fontsize=14)

fig.tight_layout()
plt.show()
fig.savefig('Growth rate simulation case2_1.png')
