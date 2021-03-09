#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 2019-09-04

"""Step_37_strain_simulate.py
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
name_list = ['100_23', '20_2', '3c6', 'ATCC53608', 'CF48_3A',
             'DSM200016', 'I5007', 'IRT', 'JCM1112',
             'LTH2584', 'LTH5448', 'MM2_3', 'MM4_1A', 'SD2112',
             'TD1', 'TMW1_112', 'TMW1_656', 'lpuph', 'mlc3', 'LR1',
             'LR10', 'LR11', 'LR12', 'LR13', 'LR14',
             'LR17', 'LR18', 'LR19', 'LR2', 'LR3',
             'LR4', 'LR6', 'LR7', 'LR8', 'LR9'
             ]

# models_df = pd.DataFrame(columns=['model_id','reaset','metset','genset'])
# for sp_name in name_list:
#
#     Lreu_tp_model = cobra.io.load_json_model(sp_name + 'Lreu_draft_3_refined_19-08-21.json')
#     Lreu_tp_model.solver = 'glpk'
#
#     keep_move_dic = {
#                     "MCOATA":"MACPMT",
#                 "DALTAL_LRE":"DALTAL",
#                 "KAS15":"kaasIII",
#                 "ALATRS":"ALATRS_1",
#                 "ASPTRS":"ASPTRS_1",
#                 "CYSTRS":"CYSTRS_1",
#                 "HISTRS":"HISTRS_1",
#                 "LEUTRS":"LEUTRS_1",
#                 "METTRS":"METTRS_1",
#                 "PHETRS":"PHETRS_1",
#                 "ARGTRS":"ARGTRS_1",
#                 "SERTRS":"SERTRS_1",
#                 "GLYTRS":"GLYTRS_1",
#                 "ILETRS":"ILETRS_1",
#                 "LYSTRS":"LYSTRS_1",
#                 "THRTRS":"THRTRS_1",
#                 "TRPTRS":"TRPTRS_1",
#                 "TYRTRS":"TYRTRS_1",
#                 "VALTRS":"VALTRS_1",
#                 "ASNTRS":"ASNTRS_1",
#
#                 "EAR40x":"BTMAT1",
#                 "PRFGS":"PRFGS_1",
#                 "PPA2":"EPPP",
#
#                 "LEUt2r":"LEUt2",
#
#                 "ARBt2":"ARAB_Lt",
#                 "ARBt2r":"ARBt2",
#                 "SERt3":"SERt2r",
#                 "SERt2r":"SERt3",
#                 #"THRt2r":"THRt3",
#
#
#                 "AGAT_LRE":"AGAT_LPL",
#                 "CLPNS_LRE":"CLPNS_LPL",
#                 "CPSS_LRE":"CPSS_LPL",
#                 #"CPSS_LRE":"CPSS_LPL2",
#                 "DAGGT_LRE":"DAGGT_LPL",
#                 "DALTAL_LPL":"DALTAL",
#                 "GLTAL_LPL":"GLTAL",
#                 "DASYN_LRE":"LPGS_LPL",
#                 "PGSA_LRE":"PGSA_LPL"
#                 }
#
#
#     for k,v in keep_move_dic.items():
#         Lreu_tp_model = My_def.merge_model.merge_reactionsid(Lreu_tp_model, k, v)
#
#      <get model information: table>
#
#     reaset = set([i.id for i in Lreu_tp_model.reactions])
#     metset = set([i.id for i in Lreu_tp_model.metabolites])
#     genset = set([i.id for i in Lreu_tp_model.genes])
#
#     nomissing_genset = set([i for i in genset if 'missing' not in i ])
#     gpr_reaset = set([i.id for i in Lreu_tp_model.reactions if i.gene_reaction_rule != ''])
#     exchenge_reaset = set([i.id for i in Lreu_tp_model.reactions if 'EX_' in i.id])
#     #gap_reaset = set([i.id for i in Lreu_draft_3_refined.reactions if 'gap' in i.notes['from']])
#     ex_metset = set([i for i in metset if '_e' in i])
#
#     # print('genes number\t',len(nomissing_genset))
#     # print('exchange\t',len(exchenge_reaset))
#     # print('inchange\t',len(reaset)-len(exchenge_reaset)-23)
#     # print('gap\t',23)
#     # print('inmetabolites\t',len(metset) -len(ex_metset) )
#     # print('exmetabolites\t',len(ex_metset))
#     #
#     locals()[sp_name+'_reaset'] = reaset - exchenge_reaset
#     locals()[sp_name+'_metset'] = metset - ex_metset
#     locals()[sp_name+'_genset'] = nomissing_genset
#
#     models_df = models_df.append({'model_id':sp_name,
#                                   'reaset':reaset - exchenge_reaset,
#                                   'metset':metset - ex_metset,
#                                   'genset':nomissing_genset}, ignore_index=True)
#
#     models_df.to_csv('models_df.tsv',sep = '\t')
#


# ax.set_xlabel('Iterations', fontname="Arial", fontsize=9)
#    ax.set_ylabel('MSE', fontname="Arial", fontsize=9)
# plt.legend(prop={'family': 'Arial', 'size': 8})
# plt.savefig('../figures/fig_b_smc_abc_fitting_curve_0203.pdf', bbox_inches='tight')

models_df_ = pd.read_csv('models_df.tsv', sep='\t', index_col=0)

for i in models_df_.index:
    sp_name = models_df_.iloc[i]['model_id']

    locals()[sp_name + '_reaset'] = eval(models_df_.iloc[i]['reaset'])
    locals()[sp_name + '_metset'] = eval(models_df_.iloc[i]['metset'])
    locals()[sp_name + '_genset'] = eval(models_df_.iloc[i]['genset'])

# %%
core_rea = locals()[name_list[0] + '_reaset']
pan_rea = locals()[name_list[0] + '_reaset']
core_met = locals()[name_list[0] + '_metset']
pan_met = locals()[name_list[0] + '_metset']

rea_count = []
met_count = []
gen_count = []

for sp_name in name_list:
    reaset = locals()[sp_name + '_reaset']
    metset = locals()[sp_name + '_metset']
    geneset = locals()[sp_name + '_genset']

    core_rea = core_rea & reaset
    pan_rea = pan_rea | reaset
    core_met = core_met & metset
    pan_met = pan_met | metset

    rea_count.append(len(reaset))
    met_count.append(len(metset))
    gen_count.append(len(geneset))

seq_gen_count = [2065, 2054, 2054, 1918, 1936, 1917, 1910, 1890, 1943, 1796, 1713, 1865,
                 1993, 2132, 1890, 1754, 1632, 1925, 1802, 2221, 2467, 2350, 2002, 2342,
                 2212, 2259, 2646, 2163, 2347, 2295, 2129, 1979, 1981, 2207, 2252]
np.array(seq_gen_count).mean()
np.array(seq_gen_count).std()

her_name = [
     'LR1',
    'LR10', 'LR11', 'LR12', 'LR13', 'LR14',
    'LR17', 'LR18', 'LR19', 'LR2', 'LR3',
    'LR4', 'LR6', 'LR7', 'LR8', 'LR9'
]

omn_name = ['JCM1112', 'MM2_3', 'MM4_1A', 'CF48_3A', 'SD2112', 'I5007', 'ATCC53608', 'DSM200016', 'IRT', 'TD1', 'mlc3',
            '100_23', '20_2', '3c6', 'lpuph',]
sour_name = ['LTH5448', 'TMW1_112', 'TMW1_656', 'LTH2584']

import pandas as pd

Lreu_d = {'sp_name': name_list, 'seq_gen_count': seq_gen_count, 'rea_count': rea_count, 'met_count': met_count,
          'gen_count': gen_count}
Lreu_df = pd.DataFrame(data=Lreu_d)
Lreu_df['group'] = 'her'

Lreu_df.loc[Lreu_df['sp_name'].isin(omn_name), ['group']] = 'omn'
Lreu_df.loc[Lreu_df['sp_name'].isin(sour_name), ['group']] = 'sou'
Lreu_df = Lreu_df.sort_values(by=['group', 'seq_gen_count'])

import matplotlib.pyplot as plt
import numpy as np

# %%
labels = Lreu_df.sp_name

# seaborn.reset_orig()
x = np.arange(len(labels))  # the label locations

# width = 0.15  # the width of the bars
plt.rcParams['figure.figsize'] = (8, 4.0)
fig, ax = plt.subplots()
plt.xlim((-0.5, 34.5))
plt.ylim((0, 3000))

rects1 = ax.plot(x, Lreu_df.seq_gen_count, '*')
rects2 = ax.plot(x, Lreu_df.rea_count, '*')
rects3 = ax.plot(x, Lreu_df.met_count, '*')
rects4 = ax.plot(x, Lreu_df.gen_count, '*')
rects5 = ax.bar((17 - 1) / 2, 4000, width=17, alpha=0.2)
rects6 = ax.bar((31 + 17 - 1) / 2, 4000, width=14, alpha=0.2)
rects7 = ax.bar((31 + 35 - 1) / 2, 4000, width=4, alpha=0.2)

plt.xticks(x, labels, rotation=90)

plt.legend(['Gens in sequence', 'Reactions', 'Metabolites', 'Gens in model', 'Herbivore', 'Omnivore', 'Sourdough', ],
           loc='lower left', mode="expand", bbox_to_anchor=(0., 1.02, 1., .102), borderaxespad=0., ncol=4)

Lreu_df[Lreu_df['group'] == 'her'].rea_count.mean()
Lreu_df[Lreu_df['group'] == 'her'].rea_count.std()

Lreu_df[Lreu_df['group'] == 'her'].met_count.mean()
Lreu_df[Lreu_df['group'] == 'her'].met_count.std()

Lreu_df[Lreu_df['group'] == 'her'].gen_count.mean()
Lreu_df[Lreu_df['group'] == 'her'].gen_count.std()

Lreu_df[Lreu_df['group'] == 'omn'].rea_count.mean()
Lreu_df[Lreu_df['group'] == 'omn'].rea_count.std()

Lreu_df[Lreu_df['group'] == 'omn'].met_count.mean()
Lreu_df[Lreu_df['group'] == 'omn'].met_count.std()

Lreu_df[Lreu_df['group'] == 'omn'].gen_count.mean()
Lreu_df[Lreu_df['group'] == 'omn'].gen_count.std()

Lreu_df[Lreu_df['group'] == 'sou'].rea_count.mean()
Lreu_df[Lreu_df['group'] == 'sou'].rea_count.std()

Lreu_df[Lreu_df['group'] == 'sou'].met_count.mean()
Lreu_df[Lreu_df['group'] == 'sou'].met_count.std()

Lreu_df[Lreu_df['group'] == 'sou'].gen_count.mean()
Lreu_df[Lreu_df['group'] == 'sou'].gen_count.std()

#
# loc='center left', bbox_to_anchor=(0.2, 1.12),ncol=3
# ax.legend(loc=2, bbox_to_anchor=(1.05,1.0),borderaxespad = 0.)     ##设置ax4中legend的位置，将其放在图外


# ax.legend(rects1, ('Gens in sequence',),rects2, ('Reactions',))
# ax.legend(rects2, ('Reactions',))
# ax.legend(rects3, ('Metabolites',))
# ax.legend(rects4, ('Gens in model',))
# ax.legend(rects5, ('Herbivore',))
# ax.legend(rects6, ('Omnivore',))
# ax.legend(rects7, ('Sourdough',))
ax.set_ylabel('Counts', )
ax.set_xlabel('GEMs features', )

fig.tight_layout()
# plt.savefig('Growth rate simulation.png')

plt.show()

# %%


colors = plt.cm.get_cmap('tab10').colors
labels = Lreu_df.sp_name

x = np.arange(len(labels))  # the label locations

# width = 0.15  # the width of the bars
plt.rcParams['figure.figsize'] = (4, 8.0)
fig, ax = plt.subplots()
plt.xlim((0, 1000))
# plt.ylim((0, 3000))

# rects1 = ax.plot(  Lreu_df.seq_gen_count,x,'*')
rects2, = ax.plot(Lreu_df.rea_count, x, '*')
rects3, = ax.plot(Lreu_df.met_count, x, '*')
rects4, = ax.plot(Lreu_df.gen_count, x, '*')
# plt.hold(True)

rects5 = ax.barh((31 + 35 - 1) / 2, 1000, height=4, alpha=0.2)
rects6 = ax.barh((31 + 16 - 1) / 2, 1000, height=15, alpha=0.2)
rects7 = ax.barh((16 - 1) / 2, 1000, height=16, alpha=0.2)

plt.yticks(x, labels, )

# legend1 = pyplot.legend(plot_lines[0], ["algo1", "algo2", "algo3"], loc=1)
# pyplot.legend([l[0] for l in plot_lines], parameters, loc=4)
# pyplot.gca().add_artist(legend1)
# legend1 = plt.legend([rects2, rects3, rects4], ['Reactions', 'Metabolites', 'Gens in model'],ncol=1,
#                     )
legend1 = plt.legend([rects2, rects3, rects4], ['Reactions', 'Metabolites', 'Gens in model'], ncol=1,
                     loc='lower right', bbox_to_anchor=(1, -0.15))
legend2 = plt.legend([rects5, rects6, rects7], ['Herbivore', 'Omnivore', 'Sourdough'], ncol=1,
                     loc='lower right', bbox_to_anchor=(0.2, -0.15))

ax.add_artist(legend1)
# plt.gca().add_artist(legend1)

# #
# loc='center left', bbox_to_anchor=(0.2, 1.12),ncol=3
# ax.legend(loc=2, bbox_to_anchor=(1.05,1.0),borderaxespad = 0.)     ##设置ax4中legend的位置，将其放在图外
# fig.savefig('fig4a1.pdf', bbox_inches='tight')

# ax.set_xlabel('Counts', )
# ax.xaxis.tick_top()
# ax.set_ylabel( 'GEMs features',)


# fig.tight_layout()
# plt.savefig('Growth rate simulation2a.png')
# plt.show()

plt.rcParams['figure.figsize'] = (3, 8.0)
fig, ax = plt.subplots()
plt.xlim((1600, 2700))
# plt.ylim((0, 3000))

rects1 = ax.plot(Lreu_df.seq_gen_count, x, '*',color = colors[5] )
rects5 = ax.barh((31 + 35 - 1) / 2, 2700, height=4, alpha=0.2)
rects6 = ax.barh((31 + 16 - 1) / 2, 2700, height=15, alpha=0.2)
rects7 = ax.barh((16 - 1) / 2, 2700, height=16, alpha=0.2)
# plt.yticks(x, labels, )
plt.yticks(x, ' ' * len(x), )

plt.legend(['Genes count (Genome))', ],
           loc='lower right', bbox_to_anchor=(0.9, -0.1))

# ax.set_xlabel('Counts', )
# ax.xaxis.tick_top()
# ax.set_ylabel( 'GEMs features',)


# fig.tight_layout()
# plt.savefig('Growth rate simulation2b.png')
# fig.savefig('fig4a2.pdf', bbox_inches='tight')

plt.show()

# %%
import My_def

her_corerea = set()
her_coremet = set()
omn_corerea = set()
omn_coremet = set()
sour_corerea = set()
sour_coremet = set()

her_panrea = set()
her_panmet = set()
omn_panrea = set()
omn_panmet = set()
sour_panrea = set()
sour_panmet = set()

a = her_corerea - omn_corerea - sour_corerea
b = omn_corerea - her_corerea - sour_corerea
c = sour_corerea - omn_corerea - her_corerea

for sp_name in her_name:
    if len(her_corerea) == 0:
        her_corerea = locals()[sp_name + '_reaset']
    else:
        her_corerea = her_corerea & locals()[sp_name + '_reaset']
        her_panrea = her_panrea | locals()[sp_name + '_reaset']

    if len(her_coremet) == 0:
        her_coremet = locals()[sp_name + '_metset']
    else:
        her_coremet = her_coremet & locals()[sp_name + '_metset']
        her_panmet = her_panmet | locals()[sp_name + '_metset']

for sp_name in omn_name:
    if len(omn_corerea) == 0:
        omn_corerea = locals()[sp_name + '_reaset']
    else:
        omn_corerea = omn_corerea & locals()[sp_name + '_reaset']
        omn_panrea = omn_panrea | locals()[sp_name + '_reaset']

    if len(omn_coremet) == 0:
        omn_coremet = locals()[sp_name + '_metset']
    else:
        omn_coremet = omn_coremet & locals()[sp_name + '_metset']
        omn_panmet = omn_panmet | locals()[sp_name + '_metset']

for sp_name in sour_name:
    if len(sour_corerea) == 0:
        sour_corerea = locals()[sp_name + '_reaset']
    else:
        sour_corerea = sour_corerea & locals()[sp_name + '_reaset']
        sour_panrea = sour_panrea | locals()[sp_name + '_reaset']

    if len(sour_coremet) == 0:
        sour_coremet = locals()[sp_name + '_metset']
    else:
        sour_coremet = sour_coremet & locals()[sp_name + '_metset']
        sour_panmet = sour_panmet | locals()[sp_name + '_metset']

print(len(her_corerea), len(her_corerea) / len(her_panrea))
print(len(her_coremet), len(her_coremet) / len(her_panmet))

print(len(omn_corerea), len(omn_corerea) / len(omn_panrea))
print(len(omn_coremet), len(omn_coremet) / len(omn_panmet))

print(len(sour_corerea), len(sour_corerea) / len(sour_panrea))
print(len(sour_coremet), len(sour_coremet) / len(sour_panmet))

from matplotlib_venn import venn3
from matplotlib import pyplot as plt

# venn3([her_corerea, omn_corerea, sour_corerea], ('Set1', 'Set2', 'Set3'))

plt.figure(figsize=(4, 4))
a = My_def.venn3_samesize([her_corerea, omn_corerea, sour_corerea], ('Herbivore', 'Omnivore', 'Sourdough'))
for text in a.set_labels:
    text.set_fontsize(14)
for text in a.subset_labels:
    text.set_fontsize(16)
plt.show()
plt.figure(figsize=(4, 4))
a = My_def.venn3_samesize([her_coremet, omn_coremet, sour_coremet], ('Herbivore', 'Omnivore', 'Sourdough'))
for text in a.set_labels:
    text.set_fontsize(14)
for text in a.subset_labels:
    text.set_fontsize(16)
plt.show()

plt.figure(figsize=(4, 4))
v = venn3(subsets=(1, 1, 1, 1, 1, 1, 1), set_labels=('Herbivore', 'Omnivore', 'Sourdough'))
textlist = (161, 27, 243, 110, 138, 34, 726)
listid = ['100', '010', '110', '001', '101', '011', '111']

plt.figure(figsize=(4, 4))
for i in range(7):
    v.get_label_by_id(listid[i]).set_text(textlist[i])
    v.get_label_by_id(listid[i]).set_fontsize(14)

# for text in v.subset_labels:
#    text.set_fontsize(16)
# plt.savefig('Venn_of_Lreu_seqs.png')
plt.show()

# %%
# Upset plot
from matplotlib import pyplot as plt
from upsetplot import from_memberships, plot, UpSet

# example = generate_counts()
# print(example)
#
# gb_degree = example.groupby(sum,group_keys=False)
# example = gb_degree.apply(lambda x: x.sort_index(ascending=False))
# print(example)
# # example.sort_values(ascending=False)
# plot(example,sort_by = 'degree',show_counts='%d',sort_categories_by =None )
# plt.show()
#
# example = example.sort_values(ascending=False)
#
#
# #
# # df = example.reorder_levels(totals.index.values)
# # agg = agg.reorder_levels(totals.index.values)
# # example['other'] = np.arange(8)
# plot(example,sort_by = 'cardinality',sort_categories_by =None,show_counts='%d')
# plt.show()

#
# plot(example)
# plt.suptitle('Ordered by degree')
# plt.show()
#
#
# plot(example, sort_by='cardinality')
# plt.suptitle('Ordered by cardinality')
# plt.show()
#
#
# plot(example, show_counts='%d')
# plt.suptitle('With counts shown')
# plt.show()
#
# plot(example, show_counts='%d', show_percentages=True)
# plt.suptitle('With counts and % shown')
# plt.show()
listid = ['100', '010', '001', '110', '101', '011', '111']
listid = [['Herbivore'],
          ['Omnivore'],
          ['Sourdough'],

          ['Herbivore', 'Omnivore'],
          ['Herbivore', 'Sourdough'],
          ['Omnivore', 'Sourdough'],
          ['Herbivore', 'Omnivore', 'Sourdough'],
          ]


def upsetplot_intreaction(textlist, title, listid, type='intersections', three_set=[], all_set=[]):
    if type == 'intersections':
        upset_plot_data = from_memberships(listid, data=textlist)
        # upset_plot_data.reorder_levels(listid[-1])
        print(upset_plot_data)

        if len(three_set) > 0:
            seven_frenquence = []
            x = []

            frenquence_i = []
            for i in three_set[0] - three_set[1] - three_set[2]:
                frenquence_i.append(all_set.count(i))
            seven_frenquence.extend(frenquence_i)
            x.extend([0] * len(frenquence_i))

            frenquence_i = []
            for i in three_set[1] - three_set[0] - three_set[2]:
                frenquence_i.append(all_set.count(i))
            seven_frenquence.extend(frenquence_i)
            x.extend([1] * len(frenquence_i))

            frenquence_i = []
            for i in three_set[2] - three_set[0] - three_set[1]:
                frenquence_i.append(all_set.count(i))
            seven_frenquence.extend(frenquence_i)
            x.extend([2] * len(frenquence_i))

            frenquence_i = []
            for i in three_set[0] & three_set[1] - three_set[2]:
                frenquence_i.append(all_set.count(i))
            seven_frenquence.extend(frenquence_i)
            x.extend([3] * len(frenquence_i))

            frenquence_i = []
            for i in three_set[0] & three_set[2] - three_set[1]:
                frenquence_i.append(all_set.count(i))
            seven_frenquence.extend(frenquence_i)
            x.extend([4] * len(frenquence_i))

            frenquence_i = []
            for i in three_set[1] & three_set[2] - three_set[0]:
                frenquence_i.append(all_set.count(i))
            seven_frenquence.extend(frenquence_i)
            x.extend([5] * len(frenquence_i))

            frenquence_i = []
            for i in three_set[1] & three_set[2] & three_set[0]:
                frenquence_i.append(all_set.count(i))
            seven_frenquence.extend(frenquence_i)
            x.extend([6] * len(frenquence_i))
            # axs = plot(upset_plot_data, sort_categories_by=None, show_counts='%d', facecolor='black', element_size=50, )

            upset_plot_data = pd.DataFrame(upset_plot_data)

            upset_plot_data['idx'] = range(7)
            print(upset_plot_data)

            temdf = pd.DataFrame({'idx': x, 'frenquency': seven_frenquence})
            upset_plot_data = upset_plot_data.merge(temdf, on='idx', how='outer', right_index=True)

            print(upset_plot_data)

            upset = UpSet(upset_plot_data, sort_categories_by=None, show_counts='%d', facecolor='black',
                          element_size=50)

            upset.add_catplot(value='frenquency', kind='strip', color='blue', )
            # upset.add_catplot(value='AGE', kind='strip', color='black')
            upset.plot()
        else:
            axs = plot(upset_plot_data, sort_categories_by=None, show_counts='%d', facecolor='black', element_size=10, )

    if type == 'union':
        textlist_union = np.array(textlist)
        textlist_initial = np.array(textlist)
        textlist_union[0] = sum(textlist_initial[[0, 3, 4, 6]])
        textlist_union[1] = sum(textlist_initial[[1, 3, 5, 6]])
        textlist_union[2] = sum(textlist_initial[[2, 4, 5, 6]])

        textlist_union[3] = sum(textlist_initial) - textlist_initial[2]
        textlist_union[4] = sum(textlist_initial) - textlist_initial[1]
        textlist_union[5] = sum(textlist_initial) - textlist_initial[0]
        textlist_union[6] = sum(textlist_initial)

        upset_plot_data = from_memberships(listid, data=textlist_union)
        axs = plot(upset_plot_data, sort_categories_by=None, show_counts='%d', facecolor='black', )
        axs['intersections'].set_ylabel('Union size')
        [i.set_alpha(0.1) for i in axs['intersections'].get_children()[:7]]

        # print(axs['intersections'].get_children())
        # axs['intersections'].bar(np.arange(7), textlist_initial, .5, color='black', label='intersections')
        axs['intersections'].bar(np.arange(7), textlist_union, .5, color='lightgrey', )

        textlist_core = np.zeros(7)

        textlist_core[3] = textlist_initial[6]
        textlist_core[4] = textlist_initial[6]
        textlist_core[5] = textlist_initial[6]
        textlist_core[6] = textlist_initial[6]
        axs['intersections'].bar(np.arange(7), textlist_core, .5, color='tab:blue', alpha=0.7)

        textlist_core_2 = np.zeros(7)

        textlist_core_2[3] += textlist_initial[3]
        textlist_core_2[4] += textlist_initial[4]
        textlist_core_2[5] += textlist_initial[5]

        axs['intersections'].bar(np.arange(7), textlist_core_2, .5, color='tab:blue',
                                 bottom=textlist_core, alpha=0.8, label='core')
        textlist_core = textlist_core + textlist_core_2

        textlist_core_3 = np.zeros(7)
        textlist_core_3[6] += textlist_initial[3] + textlist_initial[4] + textlist_initial[5]

        axs['intersections'].bar(np.arange(7), textlist_core_3, .5, color='tab:red',
                                 bottom=textlist_core, alpha=0.7, label='distinguish')
        textlist_core = textlist_core + textlist_core_3

        axs['intersections'].bar(np.arange(7), textlist_union - textlist_core, .5, color='tab:green',
                                 bottom=textlist_core, alpha=0.7, label='specific')
        textlist_core = textlist_core + textlist_core
        plt.legend(loc='upper left', bbox_to_anchor=(-0.8, 1))

        # textlist_core_4 = np.zeros(7)
        # textlist_core_4[0] = textlist_initial[0]
        # textlist_core_4[1] =  textlist_initial[1]
        # textlist_core_4[2] = textlist_initial[2]
        # textlist_core_4[3] = textlist_initial[3]
        # textlist_core_4[4] =  textlist_initial[6]
        # textlist_core_4[5] = textlist_initial[6]
        # textlist_core_4[6] = textlist_initial[6]

        axs['totals'].remove()
        print(axs)
        axs['shading'].remove()

        # axs['intersections'].get_children()[0].set_color('r')

        # c = np.array(['lightgrey'] * 21, dtype='O')
        # idx_1 = [0,4,8,9,10,12,14,16,17,18,19,20]
        # c[idx_1] = 'tab:blue'
        #
        # idx_2 = [0,]
        # c[idx_2] = 'r'
        # axs['matrix'].get_children()[0].set_facecolor(c)
        # axs['matrix'].get_children()[0].set_alpha(1)

    plt.suptitle(title)
    plt.show()


all_rea = []
all_met = []
for sp_name in name_list:
    all_rea = all_rea + list(locals()[sp_name + '_reaset'])
    all_met = all_met + list(locals()[sp_name + '_metset'])

# <Genes>

textlist = (161, 27, 110, 243, 138, 34, 726)
# upsetplot_intreaction(textlist, 'Genes', listid)
upsetplot_intreaction(textlist, 'Genes', listid, 'union')

# <Reactions>
textlist = (40, 10, 20, 95, 41, 20, 671)

# upsetplot_intreaction(textlist, 'Reactions', listid,three_set = [her_corerea, omn_corerea, sour_corerea],all_set =all_rea)
# upsetplot_intreaction(textlist, 'Reactions', listid,)
upsetplot_intreaction(textlist, 'Reactions', listid, 'union')

# <Metabolites>
textlist = (32, 4, 8, 33, 40, 12, 666)
# upsetplot_intreaction(textlist, 'Metabolites', listid,three_set =[her_coremet, omn_coremet, sour_coremet],all_set =all_met)
# upsetplot_intreaction(textlist, 'Metabolites', listid,)
upsetplot_intreaction(textlist, 'Metabolites', listid, 'union')

upsetplot_intreaction(textlist, 'Metabolites', listid, 'union')

# union

# %%
# frequency

# pan_rea = her_corerea | omn_corerea | sour_corerea
# pan_met = her_coremet | omn_coremet | sour_coremet
# print(len(pan_rea))
# print(len(pan_met))

pan_rea = her_panrea | omn_panrea | sour_panrea
pan_met = her_panmet | omn_panmet | sour_panmet
print(len(pan_rea))
print(len(pan_met))
all_rea_frenquence = []
all_met_frenquence = []

for rxn_i in pan_rea:
    all_rea_frenquence.append(all_rea.count(rxn_i))

for met_i in pan_met:
    all_met_frenquence.append(all_met.count(met_i))

her_rea_frenquence = []
omn_rea_frenquence = []
sour_rea_frenquence = []

for rxn_i in her_corerea:
    her_rea_frenquence.append(all_rea.count(rxn_i))

for rxn_i in omn_corerea:
    omn_rea_frenquence.append(all_rea.count(rxn_i))
for rxn_i in sour_corerea:
    sour_rea_frenquence.append(all_rea.count(rxn_i))

import seaborn as sns

fig, axs = plt.subplots(nrows=1, ncols=1, )

# plot violin plot
# axs.violinplot(all_rea_frenquence ,)
# axs.boxplot(all_rea_frenquence ,)
y = all_rea_frenquence + her_rea_frenquence + omn_rea_frenquence + sour_rea_frenquence
x = [0] * len(all_rea_frenquence) + [1] * len(her_rea_frenquence) + [2] * len(omn_rea_frenquence) + [3] * len(
    sour_rea_frenquence)

axs = sns.stripplot(x=x, y=y, jitter=0.1, color=".4", size=4)
axs = sns.violinplot(x=x, y=y, inner=None, color=".8")

# axs.boxplot(all_rea_frenquence ,positions = [1],)
# axs.violinplot([all_rea_frenquence , her_rea_frenquence , omn_rea_frenquence, sour_rea_frenquence]
#                ,positions = [0,1,2,3],showmeans=False, showmedians=False,showextrema=False)


axs.set_ylabel('Frequency', fontsize=14)
plt.setp(axs, xticks=[y for y in range(4)],
         xticklabels=['all', 'her', 'omn', 'sour'], )
axs.tick_params(axis='x', which='major', labelsize=14)

plt.show()

# %%
