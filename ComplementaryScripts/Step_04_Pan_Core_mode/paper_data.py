#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 3/8/21

"""paper_data.py
:description : script
:param : 
:returns: 
:rtype: 
"""

import os

import numpy as np
import pandas as pd

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

print('all genome size:\t%0.1f ± %.1f' %
      (np.array(seq_gen_count).mean(), np.array(seq_gen_count).std()))

print('all GEMs gen_count:\t%0.1f ± %.1f' %
      (Lreu_df['gen_count'].mean(), Lreu_df['gen_count'].std()))

print('all GEMs rea_count:\t%0.1f ± %.1f' %
      (Lreu_df['rea_count'].mean(), Lreu_df['rea_count'].std()))

print('all GEMs met_count:\t%0.1f ± %.1f\n' %
      (Lreu_df['met_count'].mean(), Lreu_df['met_count'].std()))

temp_df = Lreu_df[Lreu_df['group'] == 'omn']
print('omn GEMs gen_count:\t%0.1f ± %.1f' %
      (temp_df['gen_count'].mean(), temp_df['gen_count'].std()))

print('omn GEMs rea_count:\t%0.1f ± %.1f' %
      (temp_df['rea_count'].mean(), temp_df['rea_count'].std()))

print('omn GEMs met_count:\t%0.1f ± %.1f\n' %
      (temp_df['met_count'].mean(), temp_df['met_count'].std()))

temp_df = Lreu_df[Lreu_df['group'] == 'her']
print('her GEMs gen_count:\t%0.1f ± %.1f' %
      (temp_df['gen_count'].mean(), temp_df['gen_count'].std()))

print('her GEMs rea_count:\t%0.1f ± %.1f' %
      (temp_df['rea_count'].mean(), temp_df['rea_count'].std()))

print('her GEMs met_count:\t%0.1f ± %.1f\n' %
      (temp_df['met_count'].mean(), temp_df['met_count'].std()))

temp_df = Lreu_df[Lreu_df['group'] == 'sou']
print('sou GEMs gen_count:\t%0.1f ± %.1f' %
      (temp_df['gen_count'].mean(), temp_df['gen_count'].std()))

print('sou GEMs rea_count:\t%0.1f ± %.1f' %
      (temp_df['rea_count'].mean(), temp_df['rea_count'].std()))

print('sou GEMs met_count:\t%0.1f ± %.1f\n' %
      (temp_df['met_count'].mean(), temp_df['met_count'].std()))


''':in Step04_04_strains_simulate.py
print(len(her_corerea), len(her_corerea) / len(her_panrea))
print(len(her_coremet), len(her_coremet) / len(her_panmet))

print(len(omn_corerea), len(omn_corerea) / len(omn_panrea))
print(len(omn_coremet), len(omn_coremet) / len(omn_panmet))

print(len(sour_corerea), len(sour_corerea) / len(sour_panrea))
print(len(sour_coremet), len(sour_coremet) / len(sour_panmet))
847 0.859026369168357
771 0.9007009345794392

796 0.8073022312373225
715 0.8382180539273154

752 0.8191721132897604
726 0.9120603015075377


'''

textlist = (161, 27, 110, 243, 138, 34, 726)
print('specific genes:', sum(textlist[0:3]) / sum(textlist))
print('common genes:', sum(textlist[-1:]) / sum(textlist))

# <Reactions>
textlist = (40, 10, 20, 95, 41, 20, 671)
print('specific genes:', sum(textlist[0:3]) / sum(textlist))
print('common genes:', sum(textlist[-1:]) / sum(textlist))

# <Metabolites>
textlist = (32, 4, 8, 33, 40, 12, 666)
print('specific genes:', sum(textlist[0:3]) / sum(textlist))
print('common genes:', sum(textlist[-1:]) / sum(textlist))

''':
pan_rea = her_panrea | omn_panrea | sour_panrea
pan_met = her_panmet | omn_panmet | sour_panmet
print(len(pan_rea))
print(len(pan_met))
1010
870
'''
