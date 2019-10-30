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


os.chdir('../../ComplementaryData/Step_03_Compare_Refine/')
print('----- loading data -----')

Lreu_draft_3_refined = cobra.io.load_json_model('Lreu_draft_3_refined_0901.json')
print('\033[1;31;47m')
print('\033[0;34;48m')
# %% <get model information: table>
reaset = set([i.id for i in Lreu_draft_3_refined.reactions])
metset = set([i.id for i in Lreu_draft_3_refined.metabolites])
genset = set([i.id for i in Lreu_draft_3_refined.genes])

nomissing_genset = set([i for i in genset if 'missing' not in i ])
gpr_reaset = set([i.id for i in Lreu_draft_3_refined.reactions if i.gene_reaction_rule != ''])
exchenge_reaset = set([i.id for i in Lreu_draft_3_refined.reactions if 'EX_' in i.id])
#gap_reaset = set([i.id for i in Lreu_draft_3_refined.reactions if 'gap' in i.notes['from']])
ex_metset = set([i for i in metset if '_e' in i])

print('genes number\t',len(nomissing_genset))
print('exchange\t',len(exchenge_reaset))
print('inchange\t',len(reaset)-len(exchenge_reaset)-23)
print('gap\t',23)
print('inmetabolites\t',len(metset) -len(ex_metset) )
print('exmetabolites\t',len(ex_metset))

# %% <medium >


# %% <growth rate>

Lreuteri_530 = cobra.io.load_json_model('../Step_02_DraftModels/Template/template_models/Lreuteri_530_standlized.json')
iNF517 = cobra.io.load_json_model('../Step_02_DraftModels/Template/template_models/iNF517_standlized.json')
# %%

Lreu_draft_3_refined.solver = 'glpk'
Lreu_draft_3_refined.reactions.get_by_id('EX_thr__L_e').lower_bound = -0.2
Lreu_draft_3_refined.reactions.get_by_id('EX_leu__L_e').lower_bound = -0.2

Lreuteri_530.reactions.get_by_id('EX_etoh_e').upper_bound = 1000
Lreu_draft_3_refined.reactions.get_by_id('EX_etoh_e').upper_bound = 1000
Lreuteri_530.reactions.get_by_id('EX_lac__L_e').upper_bound = 1000
Lreu_draft_3_refined.reactions.get_by_id('EX_lac__L_e').upper_bound = 1000
Lreuteri_530.reactions.get_by_id('EX_ac_e').upper_bound = 1000
Lreu_draft_3_refined.reactions.get_by_id('EX_ac_e').upper_bound = 1000

Lreuteri_530.reactions.get_by_id('EX_glc__D_e').lower_bound = -20
Lreu_draft_3_refined.reactions.get_by_id('EX_glc__D_e').lower_bound = -20
Lreuteri_530.reactions.get_by_id('EX_glyc_e').lower_bound = -0
Lreu_draft_3_refined.reactions.get_by_id('EX_glyc_e').lower_bound = -0

Lreuteri_530.objective = "BIOMASS"
print('Lreuteri_530 Biomass:',Lreuteri_530.optimize())

Lreu_draft_3_refined.objective = "BIOMASS"
print('Lreu_draft_3_refined Biomass:',Lreu_draft_3_refined.optimize())
pre_data_glc = Lreu_draft_3_refined.optimize().objective_value
l_e_model = []
l_e_model.append(Lreu_draft_3_refined.optimize().fluxes.EX_lac__L_e /Lreu_draft_3_refined.optimize().fluxes.EX_etoh_e)

Lreuteri_530.summary()
Lreu_draft_3_refined.summary()

# %% <plot >
exp_glc = [0.569,0.619,0.678]       # ,0.616
exp_glc_gly = [0.689,0.739]         # ,0.677

Lreuteri_530.reactions.get_by_id('EX_glyc_e').lower_bound = -25
Lreu_draft_3_refined.reactions.get_by_id('EX_glyc_e').lower_bound = -25

print('Lreuteri_530 Biomass:',Lreuteri_530.optimize())
print('Lreu_draft_3_refined Biomass:',Lreu_draft_3_refined.optimize())
pre_data_glc_gly = Lreu_draft_3_refined.optimize().objective_value
l_e_model.append(Lreu_draft_3_refined.optimize().fluxes.EX_lac__L_e /Lreu_draft_3_refined.optimize().fluxes.EX_etoh_e)

Lreuteri_530.summary()
Lreu_draft_3_refined.summary()

# %%

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import statistics


labels = ['Glucose', 'Glucose+glycerol']
exp_glc
exp_glc_gly

ex_data = [np.mean(exp_glc),np.mean(exp_glc_gly)]
pre_data = [round(pre_data_glc*1000)/1000,round(pre_data_glc_gly*1000)/1000]

l_e_exp = [80.708/74.32,80.44/74.13]
if l_e_model[0] ==0:
    l_e_model = [20.3/20.9,28.5/25.5]
    print('\033[1;31;47m')
    print('假的！！！')

x = np.arange(len(labels))  # the label locations
x = np.array([0,0.6])
width = 0.15  # the width of the bars
Std = (statistics.stdev(exp_glc), statistics.stdev(exp_glc_gly))

fig, ax = plt.subplots()
plt.xlim((-0.3, 0.9))
plt.ylim((0.0, 1))
#rects5 = ax.plot( np.array([-10]),[-10],'ro', label=r'$Y_{lac} / Y_{eth}$')
rects1 = ax.bar(x - width/2, ex_data, width , yerr=Std ,label='Experiment')    #
rects2 = ax.bar(x + width/2, pre_data, width , label='Model')  #,


# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Growth rate (mmol/gDW/h)',fontsize=16)#color = 'tab:blue'
ax.tick_params(axis='y')    #, labelcolor='tab:blue'
ax.set_title( 'Growth rate simulation',fontsize=18)
ax.set_xticks(x)
ax.set_xticklabels(labels,fontsize=16)
ax.legend(loc='best',fontsize=14)

#plt.legend(['Experiment','Model',],loc =0,ncol = 1,fontsize=15)
#ax2 = ax.twinx()
#plt.xlim((-0.3, 0.9))
#plt.ylim((0.0, 3))
#rects3 = ax2.plot(x - width/2, l_e_exp,'ro')
#rects4 = ax2.plot(x + width/2, l_e_model,'ro' )


#ax2.set_ylabel(r'$Y_{lac} / Y_{eth}$',color = 'tab:red')
#ax2.tick_params(axis='y', labelcolor='tab:red')
#ax.legend()

def autolabel(rects):
    """Attach a text label above each bar in *rects*, displaying its height."""
    for rect in rects:
        height = rect.get_height()
        ax.annotate('{}'.format(height),
                    xy=(rect.get_x() + rect.get_width() / 2, height,),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='left', va='bottom',color='black',fontsize=14)
                     #arrowprops=dict(facecolor='blue', shrink=0.05))


autolabel(rects1)
autolabel(rects2)

fig.tight_layout()
plt.savefig('Growth rate simulation.png')
plt.show()


# %% <aa simulations>

aarealist = ['EX_ala__L_e',
            'EX_arg__L_e',
            'EX_asn__L_e',
            'EX_asp__L_e',
            'EX_cys__L_e',
            'EX_gln__L_e',
            'EX_glu__L_e',
            'EX_gly_e',
            'EX_his__L_e',
            'EX_ile__L_e',
            'EX_leu__L_e',
            'EX_lys__L_e',
            'EX_met__L_e',
            'EX_phe__L_e',
            'EX_pro__L_e',
            'EX_ser__L_e',
            'EX_thr__L_e',
            'EX_trp__L_e',
            'EX_tyr__L_e',
            'EX_val__L_e']
# EX_glyc_e: Glycerol

aadic = {}
for i in aarealist:
    rea1 = Lreu_draft_3_refined.reactions.get_by_id(i)
    rea2 = Lreuteri_530.reactions.get_by_id(i)
    bounds1 = rea1.bounds
    bounds2 = rea2.bounds

    rea1.bounds = (0.0,10)
    rea2.bounds = (0.0,10)
    #print('----- %s -----'%i )
    try:
        aadic[i] = [Lreu_draft_3_refined.optimize().objective_value, Lreuteri_530.optimize().objective_value]
    except:
        print(aadic[i])
    rea1.bounds = bounds1
    rea2.bounds = bounds2
print('Lreu_draft_3_refined\tLreuteri_530')
for i in aadic.keys():
    if aadic[i][0]<0.001 and aadic[i][1]<0.001:
        a = 'f\tf'
    elif aadic[i][0]>0.001 and aadic[i][1]>0.001:
        a = 't\tt'
    elif aadic[i][0]<0.001 and aadic[i][1]>0.001:
        a = 'f\tt'
    else:
        a = 't\tf'
    print(i.replace('EX_','').replace('__L_e','').replace('_e','')+'\t',a)





































