#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 2019-06-20

"""build_model_from_template.py
:description : script
:param : 
:returns: 
:rtype: 
"""
import cobra
import My_def
import re
def get_draft_from_template(tp_model1, blase_result_df,remove_missing_genes = True):
    '''
    build a model based on template model
    :param tp_model: template model
    :param blase_result_df: balst result
    :return:  draft model
    '''

    tp_model = tp_model1.copy()
    model = cobra.Model()
    model.description = 'GEM from template' + tp_model.id
    tp_gene_list = blase_result_df['qseqid'].to_list()
    my_gene_list = blase_result_df['sseqid'].to_list()

    for rea in tp_model.reactions:

        new_gpr_i, torf = My_def.merge_model.gpr2log(rea.gene_reaction_rule, tp_gene_list)

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
