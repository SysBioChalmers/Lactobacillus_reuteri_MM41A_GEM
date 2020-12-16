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

# %% <medium >
print('----- change medium -----')
iHL622.objective = "BIOMASS"

experiment_group = ['Glc_A', 'Glc_B', 'Glc_C', 'Glc_ave', 'Glc/Glyc_D', 'Glc/Glyc_E', 'Glc/Glyc_ave']
experiment_result = [0.57, 0.68, 0.62, 0.62, 0.74, 0.69, 0.72]
experiment_medium = {
    'EX_glc__D_e': [-17.65, -24.09, -22.23, -21.33, -19.27, -17.96, -18.62],
    'EX_glyc_e': [0.00, 0.00, 0.00, 0.00, -21.43, -18.94, -20.18],
    'EX_lac__L_e': [17.58, 24.08, 20.86, 20.84, 18.16, 16.97, 17.57],
    'EX_ac_e': [-0.19, 0.23, 0.44, 0.16, 9.29, 8.43, 8.86],
    'EX_etoh_e': [14.48, 20.05, 16.64, 17.06, 6.25, 6.46, 6.35],
    'EX_13ppd_e': [0.00, 0.00, 0.00, 0.00, 17.50, 15.62, 16.56],
    'EX_asp__L_e': [-0.51, -0.41, -0.39, -0.44, -0.50, -0.61, -0.55],
    'EX_glu__L_e': [-0.65, -0.24, -0.21, -0.37, -0.26, -0.38, -0.32],
    'EX_asn__L_e': [-0.31, 0.0, 0.0, -0.31, 0.00, 0.00, 0.0],
    'EX_ser__L_e': [-0.81, -0.56, -0.56, -0.64, -0.60, -0.68, -0.64],
    'EX_arg__L_e': [-2.14, -0.67, -0.58, -1.13, -0.36, -0.47, -0.41], }


predict_result = []
for i in range(0, len(experiment_result)):
    model = iHL622.copy()
    model.reactions.get_by_id('PFK').bounds = (0, 2)  # from LbReuteri model, the flux between PKP and EMP
    model.reactions.get_by_id('EX_ile__L_e').bounds = (0, 1000)
    model.reactions.get_by_id('EX_cys__L_e').bounds = (0, 1000)

    for rea in experiment_medium.keys():
        bound = experiment_medium[rea][i]
        if bound <= 0:
            model.reactions.get_by_id(rea).bounds = (bound, 0)
        elif bound >= 0:
            model.reactions.get_by_id(rea).bounds = (0, bound)

    predict_result.append(round(model.optimize().objective_value, 3))

experiment_result = experiment_result[0:3] + experiment_result[4:6]
predict_result = predict_result[0:3] + predict_result[4:6]
print('Experiment Biomass:', experiment_result)
print('iHL622 Biomass:', predict_result)

# %% <draw>
# df = pd.DataFrame({'Experiment':experiment_result,
#                    'Model':predict_result},)
# ax = df.plot(kind = 'bar',ylim = [0,1.0],
#              title = 'Growth rate simulation',
#              rot = 0,position = 0.5)
# ax.set_title('Growth rate simulation', fontsize=18)
# labels = ['Glucose', 'Glucose+Glycerol']
# # ax.set_xticks([1,3.5])
# # ax.set_xticklabels(labels, fontsize=16)
# plt.show()

import brewer2mpl

fig, ax = plt.subplots(figsize=(6, 4))
bmap = brewer2mpl.get_map('Set2', 'qualitative', 7)
colors = bmap.mpl_colors

plt.ylim((0.0, 1.0))
x = np.arange(0, 5)
x = np.array([0.5, 1, 1.5, 2.25, 2.75])
width = 0.15  # the width of the bars

rects1 = ax.bar(x - width / 2, experiment_result, width, label='Experiment', color=colors[1])  #
rects2 = ax.bar(x + width / 2, predict_result, width, label='Model', color=colors[0])  # ,

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Growth rate (mmol/gDW/h)', fontsize=16)  # color = 'tab:blue'
ax.tick_params(axis='y')  # , labelcolor='tab:blue'
ax.set_title('Growth rate simulation', fontsize=18)

labels = ['Glucose', 'Glucose+Glycerol']
ax.set_xticks([1, 2.5])
ax.set_xticklabels(labels, fontsize=16)
ax.legend(loc='best', fontsize=14)

fig.tight_layout()
plt.show()
fig.savefig('Growth rate simulation case1.png')
