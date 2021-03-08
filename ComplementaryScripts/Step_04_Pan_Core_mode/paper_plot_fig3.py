#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 1/25/21

"""paper_plot_fig3.py
:description : script
:param : 
:returns: 
:rtype: 
"""

import os

import matplotlib
import numpy as np
import pandas as pd

matplotlib.rc('font', family="Arial")

matplotlib.rcParams["font.family"] = 'Arial'  # 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Arial']
os.chdir('../../ComplementaryData/Step_04_Pan_Core_model/')
print('----- loading data -----')
output_table = 'models_df_products_2.tsv'

her_name = [
     'LR1',
    'LR10', 'LR11', 'LR12', 'LR13', 'LR14',
    'LR17', 'LR18', 'LR19', 'LR2', 'LR3',
    'LR4', 'LR6', 'LR7', 'LR8', 'LR9'
]

omn_name = ['JCM1112', 'MM2_3', 'MM4_1A', 'CF48_3A', 'SD2112', 'I5007', 'ATCC53608', 'DSM200016', 'IRT', 'TD1', 'mlc3',
            '100_23', '20_2', '3c6', 'lpuph',]
sour_name = ['LTH5448', 'TMW1_112', 'TMW1_656', 'LTH2584']

# models_df_products_1 = pd.read_csv(output_table, sep='\t', index_col=0)
# models_df_products_2 = pd.read_csv('models_df_products_.tsv', sep='\t', index_col=0)
# models_df_products = models_df_products_1.merge(models_df_products_2[['model_id','hista_c','dhap_c' ,'mthgxl_c', '12ppd__R_c']],
#                              how='left', on='model_id')
# models_df_products.to_csv(output_table, sep='\t')

models_df_products = pd.read_csv(output_table, sep='\t', index_col=0)

models_df_products['group'] = 'her'
models_df_products.loc[models_df_products['model_id'].isin(omn_name), ['group']] = 'omn'
models_df_products.loc[models_df_products['model_id'].isin(sour_name), ['group']] = 'sou'
models_df_products = models_df_products.sort_values(by=['group', ])

# Index(['model_id', 'growth', 'reaset', 'metset', 'genset', 'lac__L_c', 'ac_c',
#        'etoh_c', 'hista_c', 'fol_c', 'adeadocbl_c', 'ppoh_c', '13ppd_c',
#        'dhap_c', 'mthgxl_c', '12ppd__R_c', 'group'],
#       dtype='object')


products = ['lac__L_c', 'ac_c', 'etoh_c', 'hista_c', 'dhap_c', '12ppd__R_c', '13ppd_c', ]

results = {'her': [], 'omn': [], 'sou': []}
for group_i in results.keys():
    df_temp = models_df_products[(models_df_products['group'] == group_i)]
    len_i = df_temp.shape[0]

    for product_i in products:
        postive_i = df_temp[df_temp[product_i] > 0.1].shape[0] / len_i
        results[group_i].append(postive_i)
print(results)

products_plot_df = pd.DataFrame.from_dict(data=results, orient='index', columns=products)

# %%
import matplotlib.pyplot as plt

colors = plt.cm.get_cmap('Set2').colors
colors = np.array(colors)
w = 0.35

for product_i in products:
    fig, axs = plt.subplots(1, 3, figsize=(3 * w, w))
    fig.patch.set_alpha(0)
    axs[0].pie([products_plot_df[product_i][0], 1.0 - products_plot_df[product_i][0]], startangle=90,
               colors=colors[[0, -1]], )
    axs[1].pie([products_plot_df[product_i][1], 1.0 - products_plot_df[product_i][1]], startangle=90,
               colors=colors[[1, -1]])
    axs[2].pie([products_plot_df[product_i][2], 1.0 - products_plot_df[product_i][2]], startangle=90,
               colors=colors[[2, -1]])

    fig.subplots_adjust(wspace=-0.4)
    fig.savefig('fig3_' + product_i + '_percent.pdf', bbox_inches='tight')
    fig.show()

# %%
fig, axs = plt.subplots(4, 1, figsize=(w * 2, 3.5 * w))
fig.patch.set_alpha(0)
axs[0].pie([1, 0], startangle=180, labels=['Herbivore', ''],
           colors=colors[[0, -1]], textprops={'family': 'Arial', 'size': 8})
axs[1].pie([1, 0], startangle=180, labels=['Omnivore', ''],
           colors=colors[[1, -1]], textprops={'family': 'Arial', 'size': 8})
axs[2].pie([1, 0], startangle=180, labels=['Sourdough', ''],
           colors=colors[[2, -1]], textprops={'family': 'Arial', 'size': 8})
axs[3].pie([0, 1], startangle=180, labels=['', 'Negative'],
           colors=colors[[2, -1]], textprops={'family': 'Arial', 'size': 8})

fig.subplots_adjust(hspace=-0.1)
fig.savefig('fig3_' + 'lenged' + '_percent.pdf', bbox_inches='tight')
fig.show()

fig, axs = plt.subplots(1, 3, figsize=(w * 12,  w))
fig.patch.set_alpha(0)
axs[0].pie([0.5, 0.5], startangle=90, labels=[ '','Herbivore'],
           colors=colors[[-1,0 ]], textprops={'family': 'Arial', 'size': 8},)
axs[1].pie([0.5, 0.5], startangle=90, labels=[ '','Omnivore'],
           colors=colors[[ -1,1]], textprops={'family': 'Arial', 'size': 8})
axs[2].pie([0.5, 0.5], startangle=90, labels=['','Sourdough' ],
           colors=colors[[-1,2 ]], textprops={'family': 'Arial', 'size': 8})


# fig.subplots_adjust(hspace=-0.1)
fig.savefig('fig3_' + 'lenged_2' + '_percent.pdf', bbox_inches='tight')
fig.show()
