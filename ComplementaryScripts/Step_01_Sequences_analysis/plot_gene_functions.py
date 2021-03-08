#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""plot_gene_functions.py
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

matplotlib.rc('font', family="Arial")

matplotlib.rcParams["font.family"] = 'Arial'  # 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Arial']
prop = {'family': 'Arial', 'size': 8}

os.chdir('../../ComplementaryData/Step_01_Sequences_analysis/')
gene_functions_df = pd.read_csv('Lreuteri_biogaia_v03/gene_functions.tsv', sep='\t')

category = list(gene_functions_df['COG functional category'].fillna('no_COG'))
category_in_model = list(gene_functions_df['In_model'])

plot_data_df = pd.read_csv('../Initial_data/fun-20.tab', sep='\t', names=['Functional category ID',
                                                                          'RGB color',
                                                                          'description'])

plot_data = {}

for i1, i2 in zip(category, category_in_model):
    if i1 != 'no_COG':
        for ii in i1:
            if ii not in plot_data.keys():
                plot_data[ii] = [1, 0]
            else:
                plot_data[ii][0] += 1
            if i2:
                plot_data[ii][1] += 1
    else:
        if i1 not in plot_data.keys():
            plot_data[i1] = [1, 0]
        else:
            plot_data[i1][0] += 1
        if i2:
            plot_data[i1][1] += 1

plot_data_df_temp = pd.DataFrame.from_dict(plot_data, orient='index', columns=['Gene_counts', 'Gene_in_model_counts'])
plot_data_df_temp['Functional category ID'] = plot_data_df_temp.index

plot_data_df = plot_data_df.merge(plot_data_df_temp, how='outer', on='Functional category ID')

plot_data_df = plot_data_df.dropna(subset=['Gene_counts'])

# %%

# width = 0.35  # the width of the bars: can also be len(x) sequence
# labels = plot_data_df['Functional category ID']
# y1 = plot_data_df['Gene_counts']
# y2 = plot_data_df['Gene_in_model_counts']
#
# fig, ax = plt.subplots()
#
# ax.bar(labels, y2, width, label='Gene_in_model_counts')
# ax.bar(labels, y1 - y2, width, bottom=y2,
#        label='Gene_counts')
#
# ax.set_ylabel('counts', fontname="Arial", fontsize=10)
# # ax.set_title('Scores by group and gender')
# ax.legend(prop=prop)
#
# plt.show()

# %%
# width = 0.5  # the width of the bars: can also be len(x) sequence
# labels = plot_data_df['Functional category ID']
#
# x = np.arange(len(labels))
# y1 = plot_data_df['Gene_counts']
# y2 = plot_data_df['Gene_in_model_counts']
#
# fig, ax = plt.subplots()
#
# ax.bar(labels, y2, width, label='Gene_in_model_counts')
# # rects5 = ax.bar( -0.5, 100, width = 5,alpha = 0.2)
# rects6 = ax.bar((21 + 11) / 2 + 0.5, 100, width=8, alpha=0.2, label='metabolism')
# # rects7 = ax.bar( (23+20)/2-0.5,100, width = 3,alpha = 0.2)
# plt.xlim((-0.5, 23.5))
# plt.ylim((0, 100))
# plt.xticks(x, labels, rotation=90)
#
# ax.set_ylabel('Counts', fontname="Arial", fontsize=10)
# ax.set_xlabel('COG functional category ID', fontname="Arial", fontsize=10)
# ax.legend()
#
# plt.show()

# %%
# width = 0.6  # the width of the bars: can also be len(x) sequence
# labels = plot_data_df['Functional category ID']
# plot_data_df.loc[
#     plot_data_df['Functional category ID'].isin(labels[0:3]), ['group']] = 'Information Storage and Processing'
# plot_data_df.loc[
#     plot_data_df['Functional category ID'].isin(labels[3:13]), ['group']] = 'Cellular Process and Signaling'
# plot_data_df.loc[plot_data_df['Functional category ID'].isin(labels[13:21]), ['group']] = 'Metabolism'
# plot_data_df.loc[plot_data_df['Functional category ID'].isin(labels[21:]), ['group']] = 'Others'
#
# aa = plot_data_df_temp[plot_data_df_temp['Functional category ID'] == 'no_COG']['Gene_counts'].sum()
# bb = plot_data_df_temp['Gene_counts'].sum()
# 1 - aa / bb
#
# aa = plot_data_df[plot_data_df['group'] == 'Metabolism']['Gene_counts'].sum()
# bb = plot_data_df['Gene_counts'].sum()
# aa / bb
#
# aa = plot_data_df[plot_data_df['group'] == 'Cellular Process and Signaling']['Gene_counts'].sum()
# bb = plot_data_df['Gene_counts'].sum()
# aa / bb
#
# plot_data_df_temp['Gene_counts'].sum()
#
# plot_data_df_temp['Gene_in_model_counts'].sum()
# plot_data_df_temp = plot_data_df.loc[plot_data_df['Gene_in_model_counts'] > 0]
# fig, ax = plt.subplots()
#
# x_0 = 0
# x_all = np.array([])
# for group in plot_data_df_temp['group'].drop_duplicates():
#     lable_i = plot_data_df_temp[plot_data_df_temp['group'] == group]['Functional category ID']
#     x = np.arange(len(lable_i)) * 0.9 + x_0
#     x_0 = max(x) + 2
#     # y1 = plot_data_df['Gene_counts']
#     y2 = plot_data_df_temp[plot_data_df_temp['group'] == group]['Gene_in_model_counts']
#
#     group_bar = ax.bar(x, y2, width, label=group, alpha=0.8, )
#     # rects5 = ax.bar( -0.5, 100, width = 5,alpha = 0.2)
#     # rects6 = ax.bar((21 + 11) / 2 + 0.5, 100, width=8, alpha=0.2, label='metabolism')
#     # rects7 = ax.bar( (23+20)/2-0.5,100, width = 3,alpha = 0.2)
#
#     x_all = np.concatenate([x_all, x])
#
# # plt.ylim((0, 100))
# # plt.xlim((-0.5, 23.5))
#
# labels = list(plot_data_df_temp['Functional category ID'])
# labels[-1] = '*'
# plt.xticks(x_all, labels, )
#
# ax.set_ylabel('Counts', fontname="Arial", fontsize=10)
# ax.set_xlabel('COG functional category ID', fontname="Arial", fontsize=10)
# ax.legend()
# # plt.savefig('../figures/fig_b_smc_abc_fitting_curve_0203.pdf', bbox_inches='tight')
#
# plt.show()
# %%

import seaborn as sns

sns.set_style("white", )
sns.set_style("ticks", )

width = 0.6  # the width of the bars: can also be len(x) sequence
labels = plot_data_df['Functional category ID']
plot_data_df.loc[
    plot_data_df['Functional category ID'].isin(labels[0:3]), ['group']] = 'Information\nStorage &\nProcessing'
plot_data_df.loc[
    plot_data_df['Functional category ID'].isin(labels[3:13]), ['group']] = 'Cellular\nProcess &\nSignaling'
plot_data_df.loc[plot_data_df['Functional category ID'].isin(labels[13:21]), ['group']] = 'Metabolism\n\n'
plot_data_df.loc[plot_data_df['Functional category ID'].isin(labels[21:]), ['group']] = 'Others\n\n'

plot_data_df_temp = plot_data_df.loc[plot_data_df['Gene_in_model_counts'] > 0]

x = plt.cm.get_cmap('Set2')  # Set2

# import seaborn as sns
# colors = sns.color_palette("vlag")


# plt.set_cmap('jet')
# plt.set_cmap('Set2')
# matplotlib.rc("image",cmap="Set2")

fig, ax = plt.subplots(figsize=(3, 4))

x_0 = 0
x_all = np.array([])
i = 0

colors = plt.cm.get_cmap('Set2').colors
# color = [colors[1], colors[2], colors[0], colors[7]]

c1 = plt.cm.get_cmap('Set2').colors[2]
c2 = plt.cm.get_cmap('Set2').colors[1]
c3 = plt.cm.get_cmap('Set2').colors[0]
c4 = plt.cm.get_cmap('Set2').colors[7]
color = [c1, c2, c3, c4]

for group in plot_data_df_temp['group'].drop_duplicates():

    lable_i = plot_data_df_temp[plot_data_df_temp['group'] == group]['Functional category ID']
    x = np.arange(len(lable_i)) * 0.8 + x_0
    x_0 = max(x) + 1.5
    # y1 = plot_data_df['Gene_counts']
    y2 = plot_data_df_temp[plot_data_df_temp['group'] == group]['Gene_in_model_counts']

    pst = '\n {:.1%}'.format(y2.sum() / plot_data_df_temp['Gene_in_model_counts'].sum())

    group_bar = ax.barh(x, y2, width, label=group, alpha=1, color=color[i])
    # sns.catplot(x=x, y=y2, kind="bar", )

    plt.text(114, x.mean() - 0.2, group, fontsize=6, ha='center', va='center', rotation=270, weight='bold')
    plt.text(104, x.mean() - 0.2, pst, fontsize=8, ha='center', va='center', rotation=270, weight='bold')

    if 'Others' not in group:
        plt.plot([95, 130], [x_0 - 0.7, x_0 - 0.7], '--', color='black', alpha=.8, linewidth=0.8)
    # rects5 = ax.bar( -0.5, 100, width = 5,alpha = 0.2)
    # rects6 = ax.bar((21 + 11) / 2 + 0.5, 100, width=8, alpha=0.2, label='metabolism')
    # rects7 = ax.bar( (23+20)/2-0.5,100, width = 3,alpha = 0.2)

    x_all = np.concatenate([x_all, x])
    i += 1

plt.ylim((-1, 18))
plt.xlim((0, 124))

labels = list(plot_data_df_temp['Functional category ID'])
labels[-1] = '*'
plt.yticks(x_all, labels, fontname="Arial", fontsize=6)
plt.xticks(fontname="Arial", fontsize=8)

ax.invert_yaxis()
# ax.xaxis.set_label_position('top')
# ax.xaxis.set_ticks_position('top')

ax.set_xlabel('Counts', fontsize=8, family='Arial')
ax.set_ylabel('COG functional category ID', fontsize=8, family='Arial')
# legend_font = {'family': 'Arial', 'weight': 'normal', 'size': 6}
# prop = {'family': 'Arial', 'size': 5.5}

# ax.legend(loc = 'lower right',prop =prop )

# plt.title('Model genes COG functional category distribution', fontsize=8, family='Arial')
plt.savefig('fig1a_Model genes COG functional category distribution.pdf', bbox_inches='tight')

plt.show()

# %%
width = 0.6  # the width of the bars: can also be len(x) sequence
labels = plot_data_df['Functional category ID']
plot_data_df.loc[
    plot_data_df['Functional category ID'].isin(labels[0:3]), ['group']] = 'Information Storage\n and Processing'
plot_data_df.loc[
    plot_data_df['Functional category ID'].isin(labels[3:13]), ['group']] = 'Cellular Process\n and Signaling'
plot_data_df.loc[plot_data_df['Functional category ID'].isin(labels[13:21]), ['group']] = 'Metabolism'
plot_data_df.loc[plot_data_df['Functional category ID'].isin(labels[21:]), ['group']] = 'Others'

plot_data_df_temp = plot_data_df.loc[plot_data_df['Gene_counts'] > 0]
fig, ax = plt.subplots(figsize=(3, 4))

x_0 = 0
x_all = np.array([])
for group in plot_data_df_temp['group'].drop_duplicates():
    lable_i = plot_data_df_temp[plot_data_df_temp['group'] == group]['Functional category ID']
    x = np.arange(len(lable_i)) * 0.8 + x_0
    x_0 = max(x) + 1.5
    # y1 = plot_data_df['Gene_counts']
    y1 = plot_data_df_temp[plot_data_df_temp['group'] == group]['Gene_counts']

    group_bar = ax.barh(x, y1, width, label=group, alpha=0.8, )
    # rects5 = ax.bar( -0.5, 100, width = 5,alpha = 0.2)
    # rects6 = ax.bar((21 + 11) / 2 + 0.5, 100, width=8, alpha=0.2, label='metabolism')
    # rects7 = ax.bar( (23+20)/2-0.5,100, width = 3,alpha = 0.2)

    x_all = np.concatenate([x_all, x])
# 28.6%
plt.ylim((-1, 21))

# plt.xlim((0, 120))
labels = list(plot_data_df_temp['Functional category ID'])
labels[-1] = '*'
plt.yticks(x_all, labels, fontname="Arial", fontsize=6)
plt.xticks(fontname="Arial", fontsize=8)
ax.invert_yaxis()
# ax.xaxis.set_label_position('top')
# ax.xaxis.set_ticks_position('top')

ax.set_xlabel('Counts', fontsize=9, family='Arial')
ax.set_ylabel('COG functional category ID', fontsize=9, family='Arial')
legend_font = {'family': 'Arial', 'weight': 'normal', }
prop = {'family': 'Arial', 'size': 6}

ax.legend(loc='upper right', prop=prop)
# plt.title('Sequence genes COG functional category distribution_', fontsize=9, family='Arial')
plt.savefig('figS1a_Genome COG functional category distribution.pdf', bbox_inches='tight')

plt.show()

# %%

width = 0.6  # the width of the bars: can also be len(x) sequence
labels = plot_data_df['Functional category ID']
plot_data_df.loc[
    plot_data_df['Functional category ID'].isin(labels[0:3]), ['group']] = 'Information Storage and Processing'
plot_data_df.loc[
    plot_data_df['Functional category ID'].isin(labels[3:13]), ['group']] = 'Cellular Process and Signaling'
plot_data_df.loc[plot_data_df['Functional category ID'].isin(labels[13:21]), ['group']] = 'Metabolism'
plot_data_df.loc[plot_data_df['Functional category ID'].isin(labels[21:]), ['group']] = 'Others'

plot_data_df_temp = plot_data_df.loc[plot_data_df['Gene_counts'] > 0]
fig, ax = plt.subplots(figsize=(3, 4))

x_0 = 0
x_all = np.array([])
for group in plot_data_df_temp['group'].drop_duplicates():
    lable_i = plot_data_df_temp[plot_data_df_temp['group'] == group]['Functional category ID']
    x = np.arange(len(lable_i)) * 0.8 + x_0
    x_0 = max(x) + 1.5
    # y1 = plot_data_df['Gene_counts']
    y1 = plot_data_df_temp[plot_data_df_temp['group'] == group]['Gene_counts']
    y2 = plot_data_df_temp[plot_data_df_temp['group'] == group]['Gene_in_model_counts']

    # group_bar = ax.barh(np.concatenate((x,x)), np.concatenate((y1,y2)), width,label = group,alpha = 0.6)

    group1_bar = ax.barh(x, y1, width, label=group, alpha=0.7, )
    group2_bar = ax.barh(x, y2, width, alpha=0.4, color='black')

    # rects5 = ax.bar( -0.5, 100, width = 5,alpha = 0.2)
    # rects6 = ax.bar((21 + 11) / 2 + 0.5, 100, width=8, alpha=0.2, label='metabolism')
    # rects7 = ax.bar( (23+20)/2-0.5,100, width = 3,alpha = 0.2)

    x_all = np.concatenate([x_all, x])

plt.ylim((-1, 21))

# plt.xlim((0, 120))
labels = list(plot_data_df_temp['Functional category ID'])
labels[-1] = '*'
plt.yticks(x_all, labels, fontname="Arial", fontsize=6)
plt.xticks(fontname="Arial", fontsize=8)
ax.invert_yaxis()

ax.set_xlabel('Counts', fontsize=9, family='Arial')
ax.set_ylabel('COG functional category ID', fontsize=9, family='Arial')
ax.legend(loc=0, prop=prop)
# plt.title('Sequence and model genes COG category distribution', fontsize=14, family='Arial')

plt.show()
