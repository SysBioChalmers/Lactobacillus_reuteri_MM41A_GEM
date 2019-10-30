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
import cobra
import My_def


os.chdir('../../ComplementaryData/Step_03_Compare_Refine/')
print('----- loading data -----')
name_list = ['100_23','20_2','3c6','ATCC53608','CF48_3A',
             'DSM200016','I5007','IRT','JCM1112',
             'LTH2584','LTH5448','MM2_3','MM4_1A','SD2112',
             'TD1','TMW1_112','TMW1_656','lpuph','mlc3','LR1',
             'LR10','LR11','LR12','LR13','LR14',
             'LR17','LR18','LR19','LR2','LR3',
             'LR4','LR6','LR7','LR8','LR9'
            ]
for sp_name in name_list:

    Lreu_tp_model = cobra.io.load_json_model(sp_name + 'Lreu_draft_3_refined_19-08-21.json')

    keep_move_dic = {
                    "MCOATA":"MACPMT",
                "DALTAL_LRE":"DALTAL",
                "KAS15":"kaasIII",
                "ALATRS":"ALATRS_1",
                "ASPTRS":"ASPTRS_1",
                "CYSTRS":"CYSTRS_1",
                "HISTRS":"HISTRS_1",
                "LEUTRS":"LEUTRS_1",
                "METTRS":"METTRS_1",
                "PHETRS":"PHETRS_1",
                "ARGTRS":"ARGTRS_1",
                "SERTRS":"SERTRS_1",
                "GLYTRS":"GLYTRS_1",
                "ILETRS":"ILETRS_1",
                "LYSTRS":"LYSTRS_1",
                "THRTRS":"THRTRS_1",
                "TRPTRS":"TRPTRS_1",
                "TYRTRS":"TYRTRS_1",
                "VALTRS":"VALTRS_1",
                "ASNTRS":"ASNTRS_1",

                "EAR40x":"BTMAT1",
                "PRFGS":"PRFGS_1",
                "PPA2":"EPPP",

                "LEUt2r":"LEUt2",

                "ARBt2":"ARAB_Lt",
                "ARBt2r":"ARBt2",
                "SERt3":"SERt2r",
                "SERt2r":"SERt3",
                #"THRt2r":"THRt3",


                "AGAT_LRE":"AGAT_LPL",
                "CLPNS_LRE":"CLPNS_LPL",
                "CPSS_LRE":"CPSS_LPL",
                #"CPSS_LRE":"CPSS_LPL2",
                "DAGGT_LRE":"DAGGT_LPL",
                "DALTAL_LPL":"DALTAL",
                "GLTAL_LPL":"GLTAL",
                "DASYN_LRE":"LPGS_LPL",
                "PGSA_LRE":"PGSA_LPL"
                }


    for k,v in keep_move_dic.items():
        Lreu_tp_model = My_def.merge_model.merge_reactionsid(Lreu_tp_model, k, v)

    # %% <get model information: table>

    reaset = set([i.id for i in Lreu_tp_model.reactions])
    metset = set([i.id for i in Lreu_tp_model.metabolites])
    genset = set([i.id for i in Lreu_tp_model.genes])

    nomissing_genset = set([i for i in genset if 'missing' not in i ])
    gpr_reaset = set([i.id for i in Lreu_tp_model.reactions if i.gene_reaction_rule != ''])
    exchenge_reaset = set([i.id for i in Lreu_tp_model.reactions if 'EX_' in i.id])
    #gap_reaset = set([i.id for i in Lreu_draft_3_refined.reactions if 'gap' in i.notes['from']])
    ex_metset = set([i for i in metset if '_e' in i])

    # print('genes number\t',len(nomissing_genset))
    # print('exchange\t',len(exchenge_reaset))
    # print('inchange\t',len(reaset)-len(exchenge_reaset)-23)
    # print('gap\t',23)
    # print('inmetabolites\t',len(metset) -len(ex_metset) )
    # print('exmetabolites\t',len(ex_metset))
    #
    locals()[sp_name+'_reaset'] = reaset - exchenge_reaset
    locals()[sp_name+'_metset'] = metset - ex_metset
    locals()[sp_name+'_genset'] = nomissing_genset

# %%
core_rea = locals()[name_list[0]+'_reaset']
pan_rea = locals()[name_list[0]+'_reaset']
core_met = locals()[name_list[0]+'_metset']
pan_met = locals()[name_list[0]+'_metset']


rea_count = []
met_count = []
gen_count = []


for sp_name in name_list:
    reaset = locals()[sp_name+'_reaset']
    metset = locals()[sp_name+'_metset']
    geneset = locals()[sp_name+'_genset']

    core_rea = core_rea & reaset
    pan_rea = pan_rea | reaset
    core_met = core_met & metset
    pan_met = pan_met | metset

    rea_count.append(len(reaset))
    met_count.append(len(metset))
    gen_count.append(len(geneset))

seq_gen_count = [2065,2054,2054,1918,1936,1917,1910,1890,1943,1796,1713,1865,
                 1993,2132,1890,1754,1632,1925,1802,2221,2467,2350,2002,2342,
                 2212,2259,2646,2163,2347,2295,2129,1979,1981,2207,2252]

her_name = [
             'lpuph','LR1',
             'LR10','LR11','LR12','LR13','LR14',
             'LR17','LR18','LR19','LR2','LR3',
             'LR4','LR6','LR7','LR8','LR9'
            ]

omn_name = ['JCM1112','MM2_3','MM4_1A','CF48_3A','SD2112','I5007','ATCC53608','DSM200016','IRT','TD1','mlc3','100_23','20_2','3c6',]
sour_name = ['LTH5448','TMW1_112','TMW1_656','LTH2584']



import pandas as pd
Lreu_d = {'sp_name':name_list,'seq_gen_count':seq_gen_count,'rea_count':rea_count,'met_count':met_count,'gen_count':gen_count}
Lreu_df = pd.DataFrame(data=Lreu_d)
Lreu_df['group'] ='her'

Lreu_df.loc[Lreu_df['sp_name'].isin(omn_name),['group']] = 'omn'
Lreu_df.loc[Lreu_df['sp_name'].isin(sour_name),['group']] = 'sou'
Lreu_df = Lreu_df.sort_values(by=['group','seq_gen_count'])




import matplotlib.pyplot as plt
import numpy as np
import statistics

# %%
labels = Lreu_df.sp_name


x = np.arange(len(labels))  # the label locations

#width = 0.15  # the width of the bars
plt.rcParams['figure.figsize'] = (8, 4.0)
fig, ax = plt.subplots()
plt.xlim((-0.5, 34.5))
plt.ylim((0, 3000))

rects1 = ax.plot( x, Lreu_df.seq_gen_count,'*')
rects2 = ax.plot( x, Lreu_df.rea_count,'*')
rects3 = ax.plot( x, Lreu_df.met_count,'*')
rects4 = ax.plot( x, Lreu_df.gen_count,'*')
rects5 = ax.bar( (17-1)/2, 4000, width = 17,alpha = 0.2)
rects6 = ax.bar( (31+17-1)/2, 4000, width = 14,alpha = 0.2)
rects7 = ax.bar( (31+35-1)/2,4000, width = 4,alpha = 0.2)

plt.xticks(x, labels, rotation=90)

plt.legend(['Gens in sequence','Reactions','Metabolites','Gens in model','Herbivore','Omnivore','Sourdough',],
           loc ='lower left', mode="expand",bbox_to_anchor=(0., 1.02, 1., .102),borderaxespad = 0.,ncol = 4)
#
#loc='center left', bbox_to_anchor=(0.2, 1.12),ncol=3
#ax.legend(loc=2, bbox_to_anchor=(1.05,1.0),borderaxespad = 0.)     ##设置ax4中legend的位置，将其放在图外


#ax.legend(rects1, ('Gens in sequence',),rects2, ('Reactions',))
# ax.legend(rects2, ('Reactions',))
# ax.legend(rects3, ('Metabolites',))
# ax.legend(rects4, ('Gens in model',))
# ax.legend(rects5, ('Herbivore',))
# ax.legend(rects6, ('Omnivore',))
# ax.legend(rects7, ('Sourdough',))
ax.set_ylabel('Counts',)
ax.set_xlabel( 'GEMs features',)


fig.tight_layout()
plt.savefig('Growth rate simulation.png')
plt.show()

#%%
import My_def
her_corerea = set()
her_coremet = set()
omn_corerea = set()
omn_coremet = set()
sour_corerea = set()
sour_coremet = set()

for sp_name in her_name:
    if len(her_corerea)==0: her_corerea = locals()[sp_name+'_reaset']
    else:her_corerea = her_corerea & locals()[sp_name+'_reaset']

    if len(her_coremet)==0: her_coremet = locals()[sp_name+'_metset']
    else:her_coremet = her_coremet & locals()[sp_name+'_metset']

for sp_name in omn_name:
    if len(omn_corerea)==0: omn_corerea = locals()[sp_name+'_reaset']
    else:omn_corerea = omn_corerea & locals()[sp_name+'_reaset']

    if len(omn_coremet)==0: omn_coremet = locals()[sp_name+'_metset']
    else:omn_coremet = omn_coremet & locals()[sp_name+'_metset']

for sp_name in sour_name:
    if len(sour_corerea)==0: sour_corerea = locals()[sp_name+'_reaset']
    else:sour_corerea = sour_corerea & locals()[sp_name+'_reaset']

    if len(sour_coremet)==0: sour_coremet = locals()[sp_name+'_metset']
    else:omn_coremet = sour_coremet & locals()[sp_name+'_metset']



from matplotlib_venn import venn3, venn3_circles
from matplotlib import pyplot as plt

#venn3([her_corerea, omn_corerea, sour_corerea], ('Set1', 'Set2', 'Set3'))

plt.figure(figsize=(4,4))
a = My_def.venn3_samesize([her_corerea, omn_corerea, sour_corerea], ('Herbivore','Omnivore','Sourdough'))
for text in a.set_labels:
    text.set_fontsize(14)
for text in a.subset_labels:
    text.set_fontsize(16)
plt.show()
plt.figure(figsize=(4,4))
a = My_def.venn3_samesize([her_coremet, omn_coremet, sour_coremet], ('Herbivore','Omnivore','Sourdough'))
for text in a.set_labels:
    text.set_fontsize(14)
for text in a.subset_labels:
    text.set_fontsize(16)
plt.show()


plt.figure(figsize=(4,4))
v = venn3(subsets=(1, 1, 1, 1, 1, 1, 1), set_labels = ('Herbivore','Omnivore','Sourdough'))
textlist = (161, 27, 243, 110, 138, 34, 726)
listid = ['100', '010', '110', '001', '101', '011', '111']

plt.figure(figsize=(4,4))
for i in range(7):
    v.get_label_by_id(listid[i]).set_text(textlist[i])
    v.get_label_by_id(listid[i]).set_fontsize(14)

#for text in v.subset_labels:
#    text.set_fontsize(16)
plt.savefig('Venn_of_Lreu_seqs.png')
plt.show()






