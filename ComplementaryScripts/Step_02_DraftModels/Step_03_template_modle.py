#!/usr/bin/env python
#-*- coding: utf-8 -*-
# Created by lhao at 2019-02-26

import os

import cobra
import pandas as pd
from cobra import Model
import re

import My_def


def get_model_from_template(tp_model1, blase_result_df,remove_missing_genes = True):
    '''
    build a model based on template model
    :param tp_model: template model
    :param blase_result_df: balst result
    :return:  draft model
    '''

    tp_model = tp_model1.copy()
    model = Model()
    model.description = 'GEM for L.reuteri by template' + tp_model.id
    tp_gene_list = blase_result_df['qseqid'].to_list()
    my_gene_list = blase_result_df['sseqid'].to_list()

    for rea in tp_model.reactions:

        new_gpr_i, torf = My_def.gpr2log(rea.gene_reaction_rule, tp_gene_list)

        if not remove_missing_genes:
            if 'True' in new_gpr_i:
                torf = True

        if torf:
            rea.notes['from'] = [tp_model.id]
            # New gene_reaction_rule
            new_gene_reaction_rule = rea.gene_reaction_rule

            for ii in rea.genes:
                if ii.id in tp_gene_list:
                    k = tp_gene_list.index(ii.id)
                    new_gene_reaction_rule = new_gene_reaction_rule.replace(ii.id, my_gene_list[k])

            rea.gene_reaction_rule = new_gene_reaction_rule
            model.add_reactions([rea])

    removegeneslist = []
    for i in model.genes:
        if i.id not in my_gene_list:
            #i.id=i.id+'_missing'
            removegeneslist.append(i.id)

    if remove_missing_genes:
        cobra.manipulation.remove_genes(model, removegeneslist)

    else:
        for i in removegeneslist:
            reas = model.genes.get_by_id(i).reactions
            for rea in reas:
                rea.gene_reaction_rule = re.sub(i+'(?!_missing)',i+'_missing',rea.gene_reaction_rule)
        list2 = set()
        for gene in model.genes:
            if len(gene.reactions) == 0:
                list2.add(gene.id)

        cobra.manipulation.remove_genes(model, list2)

    return model


if __name__ == '__main__':
    os.chdir('../../ComplementaryData/Step_02_DraftModels/Template/')
    t_ids = ['iBT721','iNF517','iML1515']    # ['iBT721','iNF517','iMP429','iYO844','iML1515']
    for index in range(len(t_ids)):
        blast_result_df = pd.read_csv('blast/' + t_ids[index] + '_and_Lreu.csv')
        # tp_model = cobra.io.read_sbml_model('template_models/' + t_ids[index] + '.xml')
        tp_model = cobra.io.load_json_model('template_models/' + t_ids[index] + '_standlized.json')
        # Lreu_from_template_i = getmodel_from_template(tp_model, blast_result_df)BT
        Lreu_py_tp = get_model_from_template(tp_model, blast_result_df,remove_missing_genes = True)
        locals()['Lreu_from_' + t_ids[index]] = Lreu_py_tp

        cobra.io.write_sbml_model(Lreu_py_tp, 'Lreu_from_' + t_ids[index] + '.xml')
        cobra.io.save_json_model(Lreu_py_tp, 'Lreu_from_' + t_ids[index] + '.json')
        My_def.io_outtxt(Lreu_py_tp, 'Lreu_from_' + t_ids[index] + '.txt', True)

        print('done')
