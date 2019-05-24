#!/usr/bin/env python
#-*- coding: utf-8 -*-
# Created by lhao at 2019-04-24


import os
import cobra
import re
import My_def
import pandas as pd
import pickle


os.chdir('/Users/lhao/Box Sync/Project_Lreuteri/Step_03_Compare_Refine/')

with open('data.pickle', 'rb') as file:
    modelpdlist =pickle.load(file)

modelsnamelist = ['Lreu_ca.json','Lreu_ca_gp.json',
             'Lreu_se_rast_cobra.json','Lreu_se_pro_cobra.json','Lreu_se_dna_cobra.json',
             'Lreu_ra_ke.json','Lreu_ra_me.json','Lreu_ra_te.json',
             'Lreu_py_te.json',
             ]
readmodel = True
allgeneset = set()
modellist = []

for i in range(len(modelsnamelist)):
    locals()[modelsnamelist[i].split('.')[0]+ '_notepd'] = modelpdlist[i].fillna('')
    modelpd = locals()[modelsnamelist[i].split('.')[0]+ '_notepd']

    if readmodel:
        if modelsnamelist[i].endswith(r'.xml'):
            model= cobra.io.read_sbml_model(modelsnamelist[i])
        elif modelsnamelist[i].endswith(r'.json'):
            model = cobra.io.load_json_model(modelsnamelist[i])
        model.id = modelsnamelist[i].split('.')[0]
        if i not in [2,4]:
            for rea in model.reactions:
                if rea.gene_reaction_rule not in ['','Unknown']:
                    allgeneset.add(rea.notes['metanetx_id'])

        locals()[modelsnamelist[i].split('.')[0]] = model
        modellist.append(locals()[modelsnamelist[i].split('.')[0]])

mmethods = ['CarveMe', 'RAVEN_KEEG','RAVEN_MetaCyc','SEED']

gene_n = []
gene_n.append(len(Lreu_ca_gp.genes))
gene_n.append(len(Lreu_ra_ke.genes))
gene_n.append(len(Lreu_ra_me.genes))
gene_n.append(len(Lreu_se_pro_cobra.genes))
#[296, 397, 658, 516]

meta_e_n = []
meta_e_n.append(len([i.id for i in Lreu_ca_gp.metabolites if ('_e' in i.id) or (('_e0' in i.id))]))
meta_e_n.append(0)
meta_e_n.append(0)
meta_e_n.append(len([i.id for i in Lreu_se_pro_cobra.metabolites if ('_e' in i.id) or (('_e0' in i.id))]))
# [152, 0, 0, 81]

meta_c_n = []
meta_c_n.append(len([i.id for i in Lreu_ca_gp.metabolites if ('_c' in i.id) or (('_c0' in i.id))]))
meta_c_n.append(len([i.id for i in Lreu_ra_ke.metabolites ]))
meta_c_n.append(len([i.id for i in Lreu_ra_me.metabolites ]))
meta_c_n.append(len([i.id for i in Lreu_se_rast_cobra.metabolites if ('_c' in i.id) or (('_c0' in i.id))]))
#[703, 883, 846, 1090]

rea_ex_n = []
rea_ex_n.append(len([i.id for i in Lreu_ca_gp.reactions if i.id.startswith('EX_')]))
rea_ex_n.append(0)
rea_ex_n.append(0)
rea_ex_n.append(len([i.id for i in Lreu_se_rast_cobra.reactions if i.id.startswith('EX_')]))
#[150, 0, 0, 83]

rea_tra_n = []
rea_tra_n.append(len([i.id for i in Lreu_ca_gp.reactions if ("_c" in i.reaction) and ('_e' in i.reaction)]))
rea_tra_n.append(0)
rea_tra_n.append(0)
rea_tra_n.append(len([i.id for i in Lreu_se_rast_cobra.reactions  if ("_c0" in i.reaction) and ('_e0' in i.reaction)]))
#[155, 0, 0, 91]

rea_int_n = []
rea_int_n.append(len([i.id for i in Lreu_ca_gp.reactions])-rea_ex_n[0]-rea_tra_n[0])
rea_int_n.append(len([i.id for i in Lreu_ra_ke.reactions]))
rea_int_n.append(len([i.id for i in Lreu_ra_me.reactions]))
rea_int_n.append(len([i.id for i in Lreu_se_rast_cobra.reactions]) -rea_ex_n[3]-rea_tra_n[3])
#[1039, 739, 626, 916]



import numpy as np
import matplotlib.pyplot as plt

CarveMe = [gene_n[0],meta_e_n[0]]
RAVEN_KEGG =[]
RAVEN_MetaCyc = []
SEED = []


men_means  = (20, 35, 30, 35, 27)
women_means  = (25, 32, 34, 20, 25)

ind = np.arange(len(men_means))  # the x locations for the groups
width = 0.35  # the width of the bars

fig, ax = plt.subplots()
rects1 = ax.bar(ind - width/2, men_means, width,
                color='SkyBlue', label='Men')
rects2 = ax.bar(ind + width/2, women_means, width,
                color='IndianRed', label='Women')

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Scores')
ax.set_title('Scores by group and gender')
ax.set_xticks(ind)
ax.set_xticklabels(('G1', 'G2', 'G3', 'G4', 'G5'))
ax.legend()


def autolabel(rects, xpos='center'):
    """
    Attach a text label above each bar in *rects*, displaying its height.

    *xpos* indicates which side to place the text w.r.t. the center of
    the bar. It can be one of the following {'center', 'right', 'left'}.
    """

    xpos = xpos.lower()  # normalize the case of the parameter
    ha = {'center': 'center', 'right': 'left', 'left': 'right'}
    offset = {'center': 0.5, 'right': 0.57, 'left': 0.43}  # x_txt = x + w*off

    for rect in rects:
        height = rect.get_height()
        ax.text(rect.get_x() + rect.get_width()*offset[xpos], 1.01*height,
                '{}'.format(height), ha=ha[xpos], va='bottom')


autolabel(rects1, "left")
autolabel(rects2, "right")

plt.show()





