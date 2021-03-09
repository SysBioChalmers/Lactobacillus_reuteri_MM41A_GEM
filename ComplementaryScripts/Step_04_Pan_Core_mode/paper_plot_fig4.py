#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 1/20/21

"""paper_plot_fig4.py
:description : script
:param : 
:returns: 
:rtype: 
"""

import os

import matplotlib
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from upsetplot import from_memberships, plot, UpSet

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

models_df_ = pd.read_csv('models_df.tsv', sep='\t', index_col=0)

for i in models_df_.index:
    sp_name = models_df_.iloc[i]['model_id']

    locals()[sp_name + '_reaset'] = eval(models_df_.iloc[i]['reaset'])
    locals()[sp_name + '_metset'] = eval(models_df_.iloc[i]['metset'])
    locals()[sp_name + '_genset'] = eval(models_df_.iloc[i]['genset'])

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
            '100_23', '20_2', '3c6', 'lpuph', ]
sour_name = ['LTH5448', 'TMW1_112', 'TMW1_656', 'LTH2584']

Lreu_d = {'sp_name': name_list, 'seq_gen_count': seq_gen_count, 'rea_count': rea_count, 'met_count': met_count,
          'gen_count': gen_count}
Lreu_df = pd.DataFrame(data=Lreu_d)
Lreu_df['group'] = 'her'

Lreu_df.loc[Lreu_df['sp_name'].isin(omn_name), ['group']] = 'omn'
Lreu_df.loc[Lreu_df['sp_name'].isin(sour_name), ['group']] = 'sou'
Lreu_df = Lreu_df.sort_values(by=['group', 'seq_gen_count'])

np.array(Lreu_df['gen_count']).mean()
np.array(Lreu_df['gen_count']).std()

# %%
labels = Lreu_df.sp_name

colors = plt.cm.get_cmap('tab10').colors
colors = plt.cm.get_cmap('Set2').colors
colors = list(colors)
colors[3] = 'tab:blue'
# sns.set_style("ticks",)
labels = Lreu_df.sp_name
x = np.arange(len(labels))  # the label locations

# width = 0.15  # the width of the bars
fig, axs = plt.subplots(1, 2, figsize=(5, 8), sharey=True)
ax1 = axs[0]
fill_style = 'full'
rects2, = ax1.plot(Lreu_df.rea_count, labels, '*', fillstyle=fill_style,color = colors[0] )
rects3, = ax1.plot(Lreu_df.met_count, labels, '*', fillstyle=fill_style,color = colors[1])
rects4, = ax1.plot(Lreu_df.gen_count, labels, '*', fillstyle=fill_style,color = colors[2])
rects1, = ax1.plot([-100] * len(labels), labels, '*', fillstyle=fill_style, color=colors[3])

rects5 = ax1.barh((31 + 35 - 1) / 2, 1000, height=4, alpha=0.2,)
rects6 = ax1.barh((31 + 16 - 1) / 2, 1000, height=15, alpha=0.2,)
rects7 = ax1.barh((16 - 1) / 2, 1000, height=16, alpha=0.2,)
ax1.set_xlim((400, 1000))
# ax1.set_xticklabels( family='Arial', fontsize=10, )
ax1.set_yticklabels(labels, family='Arial', fontsize=8, )


# ax1.yticks(x, labels, )
#
# legend1 = plt.legend([rects2, rects3, rects4], ['Reactions', 'Metabolites', 'Gens in model'], ncol=1,
#                      loc='lower right', bbox_to_anchor=(1, -0.15))
# legend2 = plt.legend([rects5, rects6, rects7], ['Herbivore', 'Omnivore', 'Sourdough'], ncol=1,
#                      loc='lower right', bbox_to_anchor=(0.5, -0.15))
ax1.legend([rects5, rects6, rects7, rects2, rects3, rects4, rects1],
           [ 'Sourdough','Omnivore', 'Herbivore', 'Reactions', 'Metabolites', 'Genes', 'Genome size'],
           ncol=3,  # mode="expand",
           loc='lower left', bbox_to_anchor=(0, -0.14), handletextpad=2,
           prop={'family': 'Arial', 'size': 8, })
# ax.add_artist(legend1)
# plt.gca().add_artist(legend1)

# loc='center left', bbox_to_anchor=(0.2, 1.12),ncol=3
# ax.legend(loc=2, bbox_to_anchor=(1.05,1.0),borderaxespad = 0.)     ##设置ax4中legend的位置，将其放在图外
# fig.tight_layout()

ax2 = axs[1]
rects1 = ax2.plot(Lreu_df.seq_gen_count, x, '*', color=colors[3])
rects5 = ax2.barh((31 + 35 - 1) / 2, 2700, height=4, alpha=0.2)
rects6 = ax2.barh((31 + 16 - 1) / 2, 2700, height=15, alpha=0.2)
rects7 = ax2.barh((16 - 1) / 2, 2700, height=16, alpha=0.2)
# plt.yticks(x, labels, )
# plt.yticks(x, ' ' * len(x), )
ax2.set_xlim((1600, 2700))
# ax2.set_yticklabels(labels, family='Arial', fontsize=10, )

# ax2.legend(['Genes count (Genome))', ],
#            loc='lower right', bbox_to_anchor=(0.9, -0.1),
#            prop={'family': 'Arial', 'size': 8})
#
ax2.tick_params(labelsize=8)
ax1.tick_params(labelsize=8)

fig.subplots_adjust(wspace=0.06)
fig.savefig('fig4a.pdf', bbox_inches='tight')
# fig.tight_layout()

fig.show()

# %%
# Upset plot


listid = ['100', '010', '001', '110', '101', '011', '111']
listid = [['Herbivore'],
          ['Omnivore'],
          ['Sourdough'],

          ['Herbivore', 'Omnivore'],
          ['Herbivore', 'Sourdough'],
          ['Omnivore', 'Sourdough'],
          ['Herbivore', 'Omnivore', 'Sourdough'],
          ]




c1 = 'tab:blue'
c2 = 'tab:green'
c3 = 'tab:orange'
c1 = colors[2]
c2 = colors[0]
c3 = colors[1]



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
        axs = plot(upset_plot_data, sort_categories_by=None,  facecolor='black', element_size = 28)#show_counts='%d',
        axs['intersections'].set_ylabel('Union size',fontsize=10, family='Arial')
        # axs['intersections'].set_yticklabels(family='Arial', fontsize=8, )
        plt.yticks(fontname="Arial", fontsize=8)
        [i.set_alpha(0.1) for i in axs['intersections'].get_children()[:7]]

        # print(axs['intersections'].get_children())
        # axs['intersections'].bar(np.arange(7), textlist_initial, .5, color='black', label='intersections')
        axs['intersections'].bar(np.arange(7), textlist_union, .5, color='lightgrey', )

        textlist_core = np.zeros(7)

        textlist_core[3] = textlist_initial[6]
        textlist_core[4] = textlist_initial[6]
        textlist_core[5] = textlist_initial[6]
        textlist_core[6] = textlist_initial[6]
        axs['intersections'].bar(np.arange(7), textlist_core, .5, color=c1, alpha=1)

        textlist_core_2 = np.zeros(7)

        textlist_core_2[3] += textlist_initial[3]
        textlist_core_2[4] += textlist_initial[4]
        textlist_core_2[5] += textlist_initial[5]

        axs['intersections'].bar(np.arange(7), textlist_core_2, .5, color=c1,
                                 bottom=textlist_core, alpha=1, label='Common')
        textlist_core = textlist_core + textlist_core_2

        textlist_core_3 = np.zeros(7)
        textlist_core_3[6] += textlist_initial[3] + textlist_initial[4] + textlist_initial[5]

        axs['intersections'].bar(np.arange(7), textlist_core_3, .5, color=c3,
                                 bottom=textlist_core, alpha=1, label='Dispensable')
        textlist_core = textlist_core + textlist_core_3

        axs['intersections'].bar(np.arange(7), textlist_union - textlist_core, .5, color=c2,
                                 bottom=textlist_core, alpha=1, label='Specific')

        textlist_core = textlist_core + textlist_core


        # textlist_core_4 = np.zeros(7)
        # textlist_core_4[0] = textlist_initial[0]
        # textlist_core_4[1] =  textlist_initial[1]
        # textlist_core_4[2] = textlist_initial[2]
        # textlist_core_4[3] = textlist_initial[3]
        # textlist_core_4[4] =  textlist_initial[6]
        # textlist_core_4[5] = textlist_initial[6]
        # textlist_core_4[6] = textlist_initial[6]

        axs['totals'].remove()
        # print(axs)
        axs['shading'].remove()

        if title=='Metabolites':
            plt.legend(loc='center', bbox_to_anchor=(0.3, -0.65),ncol=3, handletextpad=0.4,
                       prop={'family': 'Arial', 'size': 8})
            # axs['matrix'].xaxis.set_fontsize(50)
        else:
            axs['matrix'].remove()

        # axs['intersections'].get_children()[0].set_color('r')

        # c = np.array(['lightgrey'] * 21, dtype='O')
        # idx_1 = [0,4,8,9,10,12,14,16,17,18,19,20]
        # c[idx_1] = 'tab:blue'
        #
        # idx_2 = [0,]
        # c[idx_2] = 'r'
        # axs['matrix'].get_children()[0].set_facecolor(c)
        # axs['matrix'].get_children()[0].set_alpha(1)

    # plt.suptitle(title)
    plt.savefig('fig4b_'+title+'.pdf', bbox_inches='tight')
    plt.show()


all_rea = []
all_met = []
for sp_name in name_list:
    all_rea = all_rea + list(locals()[sp_name + '_reaset'])
    all_met = all_met + list(locals()[sp_name + '_metset'])

# <Genes>

textlist = (161, 27, 110, 243, 138, 34, 726)
print('common genes: ',textlist[-1] / sum(textlist))
print('specific genes: ',(textlist[0]+textlist[1]+textlist[2] )/ sum(textlist))

# upsetplot_intreaction(textlist, 'Genes', listid)
upsetplot_intreaction(textlist, 'Genes', listid, 'union')

# <Reactions>
# textlist = (11, 10, 36, 95, 25, 20, 671)
textlist = (40, 10, 20, 95, 41, 20, 671)

print('common reactions: ',textlist[-1] / sum(textlist))
print('specific genes: ',(textlist[0]+textlist[1]+textlist[2] )/ sum(textlist))

# upsetplot_intreaction(textlist, 'Reactions', listid,three_set = [her_corerea, omn_corerea, sour_corerea],all_set =all_rea)
# upsetplot_intreaction(textlist, 'Reactions', listid,)
upsetplot_intreaction(textlist, 'Reactions', listid, 'union')

# <Metabolites>
# textlist = (2, 0, 17, 0, 11, 45, 719)
textlist = (32, 4, 8, 33, 40, 12, 666)

print('common metabolites: ',textlist[-1] / sum(textlist))
print('specific genes: ',(textlist[0]+textlist[1]+textlist[2] )/ sum(textlist))

# upsetplot_intreaction(textlist, 'Metabolites', listid,three_set =[her_coremet, omn_coremet, sour_coremet],all_set =all_met)
# upsetplot_intreaction(textlist, 'Metabolites', listid,)
upsetplot_intreaction(textlist, 'Metabolites', listid, 'union')

# upsetplot_intreaction(textlist, 'Metabolites', listid, 'union')
