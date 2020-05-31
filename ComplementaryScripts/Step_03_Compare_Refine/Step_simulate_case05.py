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

aarealist = unessential
model_result = []
model_result_gly = []
model = iHL622.copy()

for aa in aarealist:
    if model.reactions.get_by_id(aa).lower_bound == 0:
        model.reactions.get_by_id(aa).lower_bound == -1
# model.objective = 'BIOMASS'

media = {
    'EX_ala__L_e': -2.69,
    'EX_asp__L_e': -3.16,
    'EX_cys__L_e': -0.83,
    'EX_gly_e': -2.33,
    'EX_ile__L_e': -1.60,
    'EX_lys__L_e': -2.68,
    'EX_pro__L_e': -5.86,
    'EX_ser__L_e': -3.24,
}
for k, v in media.items():
    model.reactions.get_by_id(k).lower_bound = v

model.reactions.get_by_id('EX_glc__D_e').bounds = (-5, 1000)
model.reactions.get_by_id('EX_glyc_e').bounds = (0, 1000)
solution = model.optimize()
model_result.append(round(solution.objective_value, 3))

model.reactions.get_by_id('EX_glyc_e').bounds = (-5, 1000)
solution = model.optimize()
model_result_gly.append(round(solution.objective_value, 3))

for aa in aarealist:
    model_ = model.copy()

    rea = model_.reactions.get_by_id(aa)
    rea.bounds = (0, 1000)

    model_.reactions.get_by_id('EX_glyc_e').bounds = (0, 1000)
    solution = model_.optimize()
    model_result.append(round(solution.objective_value, 3))

    model_.reactions.get_by_id('EX_glyc_e').bounds = (-5, 1000)
    solution = model_.optimize()
    model_result_gly.append(round(solution.objective_value, 3))

# %%
base_ = model_result[0]
model_result = np.array(model_result) * 100 / base_
model_result_gly = np.array(model_result_gly) * 100 / base_
experiment_u = [100, 81, 98, 85, 82, 79, 96, 94, 23]
experiment_OD = [100, 97, 99, 96, 97, 34, 99, 98, 21]

experiment_u_gly = [100, 90, 96, 89, 61, 72, 94, 97, 23]
experiment_OD_gly = [100, 102, 102, 98, 94, 36, 98, 99, 19]

experiment_u_gly = np.array(experiment_u_gly) * 0.575 / 0.453
experiment_OD_gly = np.array(experiment_OD_gly) * 1.92 / 1.21

aalist = ['None'] + [i.split('_')[1] for i in aarealist]

# df = pd.DataFrame({'experiment_u': experiment_u,
#                    # 'experiment_OD': experiment_OD,
#                    'model_result': model_result,
#                    'experiment_u_gly': experiment_u_gly,
#                    # 'experiment_OD_gly': experiment_OD_gly,
#                    'model_result_gly': model_result_gly,
#                    }, index=aalist)
#
# ax = df.plot(kind = 'bar',ylim = [0,200],
#              title = 'Growth rate simulation',
#              rot = 0,position = 0.5)
# ax.set_title('Growth rate simulation', fontsize=18)
# # labels = ['Glucose', 'Glucose+Glycerol']
# # ax.set_xticks([1,3.5])
# # ax.set_xticklabels(labels, fontsize=16)
# plt.show()


# %%
import brewer2mpl

fig, ax = plt.subplots(figsize=(6, 4))
bmap = brewer2mpl.get_map('Set2', 'qualitative', 7)
colors = bmap.mpl_colors

# plt.ylim((0.0, 1.0))
x = np.arange(0, 13.5,1.5)
width = 0.2  # the width of the bars

rects1 = ax.bar(x+0.1 - width*2, experiment_u, width, label='u(max)',
                color=colors[3])  #
rects2 = ax.bar(x+0.1 - width, model_result, width, label='Model Growth',
                color=colors[2])  # ,

rects3 = ax.bar(x-0.1 + width , experiment_u_gly, width, label='u(max) with gly',
                color=colors[1])  #
rects4 = ax.bar(x-0.1 + width*2, model_result_gly, width, label='Model Growth with gly',
                color=colors[0])  # ,
ax.plot([-1, 13.5], [100, 100], "k--",alpha = 0.5)
# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Percent of None (%)', fontsize=16)  # color = 'tab:blue'
# ax.tick_params(axis='y')  # , labelcolor='tab:blue'
ax.set_ylim((0, 160))
ax.set_xlim((-1, 13))
ax.set_title('Amino acid omitted Growth rate simulation', fontsize=18)
labels = aalist
ax.set_xticks(x)
ax.set_xticklabels(labels, fontsize=14)

ax.legend(ncol=2, loc='best', fontsize=10)

fig.tight_layout()
plt.show()
fig.savefig('Growth rate simulation case5.png')


# %% plt
# import brewer2mpl
# from matplotlib.colors import LinearSegmentedColormap
# import seaborn as sns
#
# bmap = brewer2mpl.get_map('Set2', 'qualitative', 7)
# colors = bmap.mpl_colors
# # cm = LinearSegmentedColormap.from_list(
# #     'cp',[colors[1], colors[0]], N=2)
#
# sns.set()
#
# fig, ax = plt.subplots(figsize=(4, 6))
#
# im = sns.heatmap(df, linewidths=.1, ax=ax, # cmap=cm,
#                  cbar=False)
#
# cbar = plt.colorbar(im.collections[0], #orientation='horizontal',
#                     fraction=0.046, pad=0.014, shrink=0.5, aspect=10, )
#
# # cbar.set_ticks(np.array([0.25, 0.75]))
# # cbar.set_ticklabels(('no growth','growth' ))
#
#
# plt.ylabel('Amino acid omitted')
# plt.show()
# # fig.savefig('Growth rate simulation case4a.png')
#
