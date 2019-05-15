#!/usr/bin/env python
#-*- coding: utf-8 -*-
# Created by lhao at 2019-02-26


import os
import cobra
import re
from cobra import Model, Reaction, Metabolite
from Bio.Blast import NCBIXML
import My_def

os.chdir('/Users/lhao/Box Sync/Project_Lreuteri/02_DraftModels/Template/')


t_ids = ['iBT721','iNF517','iMP429','iYO844','iML1515']
t_models = ['models/'+i+'.xml'for i in t_ids]

#read blast result and have a list.

def readblast2set(file, f2t_order=True, maxE=10 ** -30, minLen=200, minIde=40):
    # set cutoff
    #   maxE              only look at genes with E-values <= this value (opt,
    #                     default 10^-30)
    #   minLen            only look at genes with alignment length >= this
    #                     value (opt, default 200)
    # 'question about minLen !!!'
    #   minIde            only look at genes with identity >= this value
    #                     (opt, default 40 (%))
    # out put : matchset, matchset.add(from_id + '\t' + to_id)

    matchset = set()

    result_handle = open(file)
    blast_records = NCBIXML.parse(result_handle)
    for blast_record in blast_records:
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                if hsp.expect <= maxE and hsp.identities >= minIde and hsp.align_length >= minLen :
                    from_id = blast_record.query.split(' ')[0]
                    to_id = alignment.hit_id
                    if f2t_order:
                        matchset.add(from_id + '\t' + to_id)
                    else:
                        matchset.add(to_id + '\t' + from_id)
    return matchset

Lreu_py_tp = Model('Lreu_py_tp')
#bigg_unmodel = cobra.io.load_json_model('/Users/lhao/Documents/Git/BIGG/universal_model.json')

allfromfromgene_list = []

for index in range(len(t_ids)):   #range(len(t_ids))
    file1 = 'Blast/' + t_ids[index] + '_in_Lreu.xml'
    file2 = 'Blast/Lreu_in_' + t_ids[index] + '.xml'

    set_ft = readblast2set(file1,False)
    set_tf = readblast2set(file2)

    matchlist = set_ft & set_tf

    targetgene_list = [j.split('\t')[1] for j in matchlist ]
    fromgene_list = [j.split('\t')[0] for j in matchlist ]
    allfromfromgene_list = allfromfromgene_list+fromgene_list

    tp_model = cobra.io.read_sbml_model(t_models[index])

    if index in [0,2]:
        for met in tp_model.metabolites:
            met.id = met.id.replace('LSQBKT', '')
            met.id = met.id.replace('_RSQBKT', '')
    #Lreu_from_1363_cobrapy = cobra.io.read_sbml_model("Lreu_from_1363_cobrapy.xml")

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
                if i.reaction !=rea.reaction or i.bounds != rea.bounds:
                    i.id = i.id + '_' + t_ids[index]
                    print(i,i.bounds)
                    print(rea,rea.bounds)
                    Lreu_py_tp.add_reactions([i])
                else:
                    rea.notes['from'] = rea.notes['from'] + i.notes['from']

                    if rea.gene_reaction_rule !=i.gene_reaction_rule:
                        print (rea.gene_reaction_rule,i.gene_reaction_rule)
                        rea.gene_reaction_rule = '(' + rea.gene_reaction_rule + ') or (' + i.gene_reaction_rule + ')'
            else:
                Lreu_py_tp.add_reactions([i])



removegeneslist = [i for i in Lreu_py_tp.genes if i.id not in allfromfromgene_list]
cobra.manipulation.remove_genes(Lreu_py_tp,removegeneslist)

#cobra.io.write_sbml_model(Lreu_py_tp, 'Lreu_py_te.xml')
#My_def.io_outtxt(Lreu_py_tp,"Lreu_py_te.txt",True)

print ('done')


'''
Lreu_ra_tp = cobra.io.read_sbml_model('Lreu_ra_te.xml')
for i in Lreu_ra_tp.metabolites:
    if '_i' in i.id:
        for ii in i.reactions:
            ii.reaction = ii.reaction.replace(i.id,i.id.split('_i')[0])

My_def.io_outtxt(Lreu_ra_tp,"Lreu_ra_te.txt",True)

'''



