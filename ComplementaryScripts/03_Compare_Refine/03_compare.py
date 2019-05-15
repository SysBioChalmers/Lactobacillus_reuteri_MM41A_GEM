#!/usr/bin/env python
#-*- coding: utf-8 -*-
# Created by lhao at 2019-04-19

import os
import cobra
import re
import My_def
import pandas as pd
import pickle
from matplotlib import pyplot as plt
from matplotlib_venn import venn3, venn3_circles
from matplotlib_venn import venn2, venn2_circles

os.chdir('/Users/lhao/Box Sync/Project_Lreuteri/03_Compare_Refine/')

with open('data.pickle', 'rb') as file:
    modelpdlist =pickle.load(file)

modelsnamelist = ['Lreu_ca.json','Lreu_ca_gp.json',
             'Lreu_se_rast.json','Lreu_se_pro.json','Lreu_se_dna.json',
             'Lreu_ra_ke.json','Lreu_ra_me.json','Lreu_ra_te.json',
             'Lreu_py_te.json',
             ]
readmodel = True

for i in range(len(modelsnamelist)):
    locals()[modelsnamelist[i].split('.')[0]+ '_notepd'] = modelpdlist[i].fillna('')
    modelpd = locals()[modelsnamelist[i].split('.')[0]+ '_notepd']
    if readmodel:
        if modelsnamelist[i].endswith(r'.xml'):
            locals()[modelsnamelist[i].split('.')[0] ]= cobra.io.read_sbml_model(modelsnamelist[i])
        elif modelsnamelist[i].endswith(r'.json'):
            locals()[modelsnamelist[i].split('.')[0]] = cobra.io.load_json_model(modelsnamelist[i])

    if False:

        modelpd['rea_set'][modelpd['rea_set'].str.startswith('EX_')] = ''
        modelpd['rea_set'][modelpd['rea_set'].str.startswith('sink_')] = ''

        modelpd['gene_set'].to_csv(modelsnamelist[i].split('.')[0]+'_gene.csv',index = False)
        modelpd['rea_set'].to_csv(modelsnamelist[i].split('.')[0] + '_rea.csv', index=False)
        modelpd['met_set'].to_csv(modelsnamelist[i].split('.')[0] + '_met.csv', index=False)


if False:
    # seed compare
    figure, axes = plt.subplots(1, 3)
    axes[0].set_title("gene")
    axes[1].set_title("rea")
    axes[2].set_title("met")

    fg1 = My_def.venn3_samesize([set(Lreu_se_pro_notepd.gene_set),
                 set(Lreu_se_dna_notepd.gene_set),
                 set(Lreu_se_rast_notepd.gene_set)],
                ('pro', 'dna','rast'), ax=axes[0])

    fg2 = My_def.venn3_samesize([set(Lreu_se_pro_notepd.rea_set),
                 set(Lreu_se_dna_notepd.rea_set),
                 set(Lreu_se_rast_notepd.rea_set)],
                ('pro', 'dna','rast'), ax=axes[1])

    fg3 = My_def.venn3_samesize([set(Lreu_se_pro_notepd.met_set),
                 set(Lreu_se_dna_notepd.met_set),
                 set(Lreu_se_rast_notepd.met_set)],
                ('pro', 'dna','rast'), ax=axes[2])
    plt.show()

if False:
    ##CarveMe compare
    figure, axes = plt.subplots(1, 3)
    axes[0].set_title("gene")
    axes[1].set_title("rea")
    axes[2].set_title("met")
    fg1= venn2([set(Lreu_ra_ke_notepd.gene_set),set(Lreu_ra_me_notepd.gene_set)],
               ('ke','me'),ax=axes[0])


    #fg1.get_patch_by_id('10').set_color('Aquamarine')

    fg2= venn2([set(Lreu_ra_ke_notepd.rea_set),set(Lreu_ra_me_notepd.rea_set)],
               ('ke','me'),ax=axes[1])


    fg3= venn2([set(Lreu_ra_ke_notepd.met_set),set(Lreu_ra_me_notepd.met_set)],
               ('ke','me'),ax=axes[2])

    plt.show()


if False:
    ##CarveMe compare
    figure, axes = plt.subplots(1, 3)
    axes[0].set_title("gene")
    axes[1].set_title("rea")
    axes[2].set_title("met")
    fg1= venn2([set(Lreu_ca_notepd.gene_set),set(Lreu_ca_gp_notepd.gene_set)],
               ('ca','ca_gp'),ax=axes[0])
    fg1.get_patch_by_id('11').set_alpha(0.5)
    fg1.get_patch_by_id('11').set_color('#80BF40')

    #fg1.get_patch_by_id('10').set_color('Aquamarine')

    fg2= venn2([set(Lreu_ca_notepd.rea_set),set(Lreu_ca_gp_notepd.rea_set)],
               ('ca','ca_gp'),ax=axes[1])
    fg2.get_patch_by_id('11').set_alpha(0.5)
    fg2.get_patch_by_id('11').set_color('#80BF40')

    fg3= venn2([set(Lreu_ca_notepd.met_set),set(Lreu_ca_gp_notepd.met_set)],
               ('ca','ca_gp'),ax=axes[2])
    fg3.get_patch_by_id('11').set_alpha(0.5)
    fg3.get_patch_by_id('11').set_color('#80BF40')
    plt.show()

if False:
    ##template compare
    figure, axes = plt.subplots(1, 3)
    axes[0].set_title("gene")
    axes[1].set_title("rea")
    axes[2].set_title("met")
    fg1= venn2([set(Lreu_ra_te_notepd.gene_set),set(Lreu_py_te_notepd.gene_set)],
               ('ra_te','py_te'),ax=axes[0])


    #fg1.get_patch_by_id('10').set_color('Aquamarine')

    fg2= venn2([set(Lreu_ra_te_notepd.rea_set),set(Lreu_py_te_notepd.rea_set)],
               ('ra_te','py_te'),ax=axes[1])


    fg3= venn2([set(Lreu_ra_te_notepd.met_set),set(Lreu_py_te_notepd.met_set)],
               ('ra_te','py_te'),ax=axes[2])

    plt.show()



