#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 2019-05-27

"""Step_01_model_comparison.py
:description : script
:param : 
:returns: 
:rtype: 
"""
import os

import cobra
from matplotlib import pyplot as plt
from matplotlib_venn import venn2
import pandas as pd
import My_def
from My_def.model_report import *


if  __name__ == '__main__':
    os.chdir('../../ComplementaryData/Step_02_DraftModels/')
    # %% <load data>
    Lreu_ca = cobra.io.load_json_model('CarveMe/Lreu_ca.json')
    Lreu_ca_gp = cobra.io.load_json_model('CarveMe/Lreu_ca_gp.json')
    Lreu_from_iNF517 = cobra.io.load_json_model('Template/Lreu_from_iNF517.json')
    Lreu_from_iBT721 = cobra.io.load_json_model('Template/Lreu_from_iBT721.json')
    Lreu_from_iML1515 = cobra.io.load_json_model('Template/Lreu_from_iML1515.json')
    bigg_rea_df = pd.read_csv('../bigg_database/bigg_rea_df.csv', sep='\t')
    bigg_met_df = pd.read_csv('../bigg_database/bigg_met_df.csv', sep='\t')

    Lreu_ca_genes = [i.id for i in Lreu_ca.genes]
    Lreu_ca_gp_genes = [i.id for i in Lreu_ca_gp.genes]
    Lreu_ca_reas = [i.id for i in Lreu_ca.reactions]
    Lreu_ca_gp_reas = [i.id for i in Lreu_ca_gp.reactions]
    Lreu_ca_mets = [i.id for i in Lreu_ca.metabolites]
    Lreu_ca_gp_mets = [i.id for i in Lreu_ca_gp.metabolites]

    # %% <fig compare Lreu_ca and Lreu_ca_gp>
    # Lreu_ca_gp have more
    figure, axes = plt.subplots(1, 3)
    axes[0].set_title("gene")
    axes[1].set_title("rea")
    axes[2].set_title("met")

    fg1 = venn2([set(Lreu_ca_genes), set(Lreu_ca_gp_genes)],
                ('Normal','Gram positive' ), ax=axes[0])

    # fg1.get_patch_by_id('10').set_color('Aquamarine')

    fg2 = venn2([set(Lreu_ca_reas), set(Lreu_ca_reas)],
                ('Normal','Gram positive'), ax=axes[1])

    fg3 = venn2([set(Lreu_ca_mets), set(Lreu_ca_gp_mets)],
                ('Normal','Gram positive'), ax=axes[2])
    plt.show()

    Lreu_from_iBT721_genes = [i.id for i in Lreu_from_iBT721.genes]
    Lreu_from_iBT721_reas = [i.id for i in Lreu_from_iBT721.reactions]
    Lreu_from_iBT721_mets = [i.id for i in Lreu_from_iBT721.metabolites]
    Lreu_from_iNF517_genes = [i.id for i in Lreu_from_iNF517.genes]
    Lreu_from_iNF517_reas = [i.id for i in Lreu_from_iNF517.reactions]
    Lreu_from_iNF517_mets = [i.id for i in Lreu_from_iNF517.metabolites]
    Lreu_from_iML1515_genes = [i.id for i in Lreu_from_iML1515.genes]
    Lreu_from_iML1515_reas = [i.id for i in Lreu_from_iML1515.reactions]
    Lreu_from_iML1515_mets = [i.id for i in Lreu_from_iML1515.metabolites]

    # %% <fig compare templated based method models and Lreu_ca_gp>
    # just a overview
    figure_2, axes = plt.subplots(1, 3)
    axes[0].set_title("gene")
    axes[1].set_title("rea")
    axes[2].set_title("met")

    fg_1 = My_def.venn3_samesize([set(Lreu_from_iBT721_genes),
                 set(Lreu_from_iNF517_genes),
                 set(Lreu_from_iML1515_genes)],
                ('iBT721', 'iNF517','iML1515'), ax=axes[0])

    fg_2 = My_def.venn3_samesize([set(Lreu_from_iBT721_reas),
                 set(Lreu_from_iNF517_reas),
                 set(Lreu_from_iML1515_reas)],
                ('iBT721', 'iNF517','iML1515'), ax=axes[1])

    fg_3 = My_def.venn3_samesize([set(Lreu_from_iBT721_mets),
                 set(Lreu_from_iNF517_mets),
                 set(Lreu_from_iML1515_mets)],
                ('iBT721', 'iNF517','iML1515'), ax=axes[2])
    plt.show()




