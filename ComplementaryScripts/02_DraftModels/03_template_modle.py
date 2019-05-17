#!/usr/bin/env python
#-*- coding: utf-8 -*-
# Created by lhao at 2019-02-26


import os
import cobra
import re
import pandas as pd
from cobra import Model, Reaction, Metabolite
import My_def

os.chdir('../../ComplementaryData/02_DraftModels/Template/')


t_ids = ['iBT721','iNF517']     #['iBT721','iNF517','iMP429','iYO844','iML1515']


#read blast result and have a list.

Lreu_py_tp = Model('Lreu_py_tp')
Lreu_py_tp.description = 'GEM for L.reuteri by "iBT721","iNF517" template'
#bigg_unmodel = cobra.io.load_json_model('/Users/lhao/Documents/Git/BIGG/universal_model.json')

allfromfromgene_list = []

for index in range(len(t_ids)):   #range(len(t_ids))

    result_df = pd.read_csv('blast/' + t_ids[index] + '_and_Lreu.csv')

    matchlist = result_df['qseqid']

    targetgene_list = result_df['qseqid'].to_list()
    fromgene_list = result_df['sseqid'].to_list()
    allfromfromgene_list = allfromfromgene_list+fromgene_list

    tp_model = cobra.io.read_sbml_model('template_models/' + t_ids[index] + '.xml')

    # standardize iBT721 and iMP429 mets id
    if index in [0,2]:
        for met in tp_model.metabolites:
            met.id = met.id.replace('LSQBKT', '')
            met.id = met.id.replace('_RSQBKT', '')

    for i in tp_model.reactions:
        new_gpr_i, torf = My_def.gpr2log(i.gene_reaction_rule, targetgene_list)
        i.notes['from'] = [t_ids[index]]

        if torf:

            # New gene_reaction_rule
            new_gene_reaction_rule = i.gene_reaction_rule
            for ii  in i.genes:
                if ii.id in targetgene_list:
                    k = targetgene_list.index(ii.id)
                    new_gene_reaction_rule = new_gene_reaction_rule.replace(ii.id,fromgene_list[k])
            i.gene_reaction_rule = new_gene_reaction_rule

            reaset = set([rea.id for rea in Lreu_py_tp.reactions])

            if i.id in reaset:
                rea = Lreu_py_tp.reactions.get_by_id(i.id)

                if i.reaction !=rea.reaction or i.bounds != rea.bounds  or rea.gene_reaction_rule !=i.gene_reaction_rule:
                    i.id = i.id + '_' + t_ids[index]
                    #print(i,i.bounds)
                    #print(rea,rea.bounds)
                    # print (rea.gene_reaction_rule,'********',i.gene_reaction_rule)
                    Lreu_py_tp.add_reactions([i])
                else:
                    rea.notes['from'] = rea.notes['from'] + i.notes['from']
            else:
                Lreu_py_tp.add_reactions([i])


removegeneslist = [i for i in Lreu_py_tp.genes if i.id not in allfromfromgene_list]

cobra.manipulation.remove_genes(Lreu_py_tp,removegeneslist)

cobra.io.write_sbml_model(Lreu_py_tp, '../../../ModelFiles/Lreu_py_te.xml')
cobra.io.save_json_model( Lreu_py_tp, '../../../ModelFiles/Lreu_py_te.json')
My_def.io_outtxt(Lreu_py_tp,"../../../ModelFiles/Lreu_py_te.txt",True)

print ('done')






