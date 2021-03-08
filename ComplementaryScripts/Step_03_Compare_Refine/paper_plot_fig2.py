#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 1/19/21

"""paper_plot_fig2.py
:description : script
:param : 
:returns: 
:rtype: 
"""

import os

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap

matplotlib.rc('font', family="Arial")

matplotlib.rcParams["font.family"] = 'Arial'  # 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Arial']
prop = {'family': 'Arial', 'size': 8}

os.chdir('../../ComplementaryData/Step_03_Compare_Refine/')

# %%a

experiment_result = [0.57, 0.68, 0.62, 0.74, 0.69]
predict_result = [0.653, 0.744, 0.667, 0.866, 0.815]

colors = sns.color_palette("Set2")
c1 = colors[2]
c2 = colors[0]

c1 = colors[-1]
c2 = colors[2]

c1 = plt.cm.get_cmap('Set2').colors[-1]  # grey
c2 = plt.cm.get_cmap('Set2').colors[0]
# c2 = plt.cm.get_cmap('tab20c').colors[1]
c3 = plt.cm.get_cmap('Set2').colors[1]
c4 = plt.cm.get_cmap('Set2').colors[0]
# c3 = plt.cm.get_cmap('Dark2').colors[0]
# c4 = plt.cm.get_cmap('Dark2').colors[1]

# c2 = plt.cm.get_cmap('Set2').colors[0]
# # c2 = plt.cm.get_cmap('tab20c').colors[1]
# c3 = plt.cm.get_cmap('Set2').colors[2]
# c4 = plt.cm.get_cmap('Set2').colors[1]


sns.set_style("ticks", )
fig, ax = plt.subplots(figsize=(3.5, 3))
k = 0.5
ax.spines['left'].set_linewidth(k)
ax.spines['right'].set_linewidth(k)
ax.spines['bottom'].set_linewidth(k)
ax.spines['top'].set_linewidth(k)
x = np.arange(0, 5)

x = np.array([0.5, 1, 1.5, 2.25, 2.75])
width = 0.2  # the width of the bars

rects1 = ax.bar(x - width / 2, experiment_result, width, label='Experiment',
                color=c3, )  # edgecolor='black', linewidth=1.2
rects2 = ax.bar(x + width / 2, predict_result, width, label='Model', color=c4, )  # ,

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Growth rate ($h^{-1}$)', fontsize=10, family='Arial', )  # color = 'tab:blue'
# ax.tick_params(axis='y')  # , labelcolor='tab:blue'
# ax.set_title('Growth rate simulation', fontsize=12, family='Arial')

labels = ['Glucose', 'Glucose+Glycerol']
ax.set_xticks([1, 2.5])
ax.set_xticklabels(labels, fontsize=10, family='Arial')
plt.yticks(fontname="Arial", fontsize=8)
ax.legend(loc='best', prop={'family': 'Arial', 'size': 8})
ax.tick_params(width=0.75)

plt.ylim((0.0, 1.0))
# ax.legend(loc = 'lower right',prop =prop )
fig.tight_layout()
fig.savefig('fig2a_Growth rate simulation.pdf', bbox_inches='tight')
fig.show()

# %%

cm = LinearSegmentedColormap.from_list(
    'cp', [c1, c2], N=2)
df = pd.read_csv('amino_acid.tsv', sep='\t', index_col=0)
aalist = df.index

sns.set()
fig, ax = plt.subplots(figsize=(8, 3))

im = sns.heatmap(df.T, linewidths=.5, ax=ax, cmap=cm,
                 cbar=False, square=True)

cbar = plt.colorbar(im.collections[0],  # orientation='horizontal',
                    fraction=0.046, pad=0.014, shrink=0.5, aspect=10, )

cbar.set_ticks(np.array([0.25, 0.75]))
cbar.set_ticklabels(('not\ngrowth', 'growth'), )
cbar.ax.tick_params(length=0, labelsize=10)

x = plt.xticks()[0]
plt.xticks(x, aalist, rotation=45, fontsize=10, family='Arial')
plt.yticks(fontsize=10, family='Arial')
ax.tick_params(length=0, )
fig.tight_layout()
fig.savefig('fig2b_amino_acid_simulation.pdf', bbox_inches='tight')
fig.show()

# %%
prd_df = pd.read_csv('products_df.tsv', sep='\t', index_col=0)
prd_list = prd_df.index
prd_list = ['lactate', 'acteate', 'ethanol ', 'histamine', 'folate', 'cobalamin',
            '1-propanol', '1,3-ppd']

sns.set()

fig, ax = plt.subplots(figsize=(4.5, 3))

im = sns.heatmap(prd_df.T, linewidths=.5, ax=ax, cmap=cm,
                 cbar=False, square=True)

cbar = plt.colorbar(im.collections[0],  # orientation='horizontal',
                    fraction=0.046, pad=0.014, shrink=0.5, aspect=10, )

cbar.set_ticks(np.array([0.25, 0.75]))
cbar.set_ticklabels(('not\nproduce', 'produce'))
cbar.ax.tick_params(length=0, labelsize=10)

x = plt.xticks()[0]
plt.xticks(x, prd_list, rotation=60, fontsize=10, family='Arial')
plt.yticks(fontsize=10, family='Arial')
ax.tick_params(length=0)
fig.tight_layout()
fig.savefig('fig2c_products_simulation.pdf', bbox_inches='tight')

fig.show()
