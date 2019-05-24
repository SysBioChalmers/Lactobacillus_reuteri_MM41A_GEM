#!/usr/bin/env python
#-*- coding: utf-8 -*-
# Created by lhao at 2019-04-22

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
biggmodel = cobra.io.load_json_model('/Users/lhao/Documents/Git/BIGG/universal_model.json')
if True:
    #addrea
    ca = set(Lreu_ca_gp_notepd.rea_set)
    ke = set(Lreu_ra_ke_notepd.rea_set)
    me = set(Lreu_ra_me_notepd.rea_set)
    se = set(Lreu_se_pro_cobra_notepd.rea_set)

    rt  = set(Lreu_ra_te_notepd.rea_set)
    pt  = set(Lreu_py_te_notepd.rea_set)
    iNF517 = cobra.io.load_json_model('iNF517.json')
    iNF517realist = [i.id for i in iNF517.reactions]
    iNF517mnxid = My_def.mapIDsViaMNXref('rxns', iNF517realist, 'bigg', 'metanetx')[1]

    tempmodel = My_def.merge_model(Lreu_ra_te,Lreu_py_te)
    print(len(tempmodel.reactions),len(tempmodel.genes))
    #327
    add_three = True
    add_two = True
    if add_three:
        allset = ((ca & ke & me ) | (ca & ke  & se)| (ca & me & se))-(rt|pt)
        allset2 = (ke & me & se)-(rt|pt) - allset
    if add_two:
        allset = (allset|(ca & ke)|(ca & me)|(ca & se))-(rt|pt)
        allset2 = (allset2|(ke & me)|(ke & se)|(se & me))-(rt|pt) - allset
    gaplist = []
    addlist = []

    #allbigglist = My_def.mapIDsViaMNXref('rxns',allset,'metanetx','bigg')[0]
    for i in Lreu_ca_gp.reactions:
        if i.notes['metanetx_id'] in allset:
            i.notes['from'] = ['Lreu_ca_gp']
            tempmodel.add_reaction(i);
    print(len(tempmodel.reactions),len(tempmodel.genes))
    #662

    for i in allset2&set(iNF517mnxid):
        reaindex = iNF517mnxid.index(i)

        rea = iNF517.reactions[reaindex]
        rea.gene_reaction_rule = ''
        rea.notes['metanetx_id'] = i
        rea.notes['from'] = ['iNF517']
        tempmodel.add_reaction(rea)
    allset2 = allset2-set(iNF517mnxid)

    allset2 = list(allset2)
    allset2list = My_def.mapIDsViaMNXref('rxns',allset2,'metanetx','bigg')[0]
    Repeatedrealist = []

    for i in range(len(allset2)):
        rea = ''
        if allset2list[i] ==['']:
            #print("option1",allset2[i])
            gaplist.append(allset2[i])
            continue
        elif len(allset2list[i]) ==1:
            rea = biggmodel.reactions.get_by_id(allset2list[i][0])
            #print("option2",'maybe is this rea')
        else:
            ###!!!!!!!

            biggids = []
            reaqeus = []
            for bigg_i in allset2list[i]:
                if not (bigg_i.startswith('R_')):
                    try :
                        reaqeu = biggmodel.reactions.get_by_id(bigg_i).reaction
                        if ('_c' in reaqeu) or ('_e' in reaqeu):
                            reaqeus.append(reaqeu)
                            biggids.append(bigg_i)
                    except :
                        pass
            if len(set(reaqeus)) ==1:
                rea = biggmodel.reactions.get_by_id(biggids[0])

            elif len(set(reaqeus)) == 2:

                nc = 0
                for kk in set(reaqeus):
                    if '_c' in kk:
                        nc +=1
                if nc == 1:


                #if ('_c' in list[set(reaqeus)][0] and '_c' not in list[set(reaqeus)][1]) or ('_c' in list[set(reaqeus)][1] and '_c' not in list[set(reaqeus)][0]):

                    for kk in biggids:

                        if '_c ' in biggmodel.reactions.get_by_id(kk).reaction:
                            rea = biggmodel.reactions.get_by_id(kk)
                            break

                else:
                    Repeatedrealist.append(biggids)
                    print('duogehouxuan!!!')
                    continue

            elif len(set(reaqeus)) >2:

                Repeatedrealist.append(biggids)

                print('duogehouxuan!!!')
                continue
            elif len(set(reaqeus)) ==0:
                print("option3", allset2[i])
                gaplist.append(allset2[i])
                continue

            #print("option3")

        if rea !='':
            rea.notes['metanetx_id'] = allset2[i]

            tempmodel.add_reaction(rea);
            rea.notes['from'] = []

            addlist.append(allset2[i])

    print(len(tempmodel.reactions),len(tempmodel.genes))



if True:
    #add gene
    ridlist = [i.id for i in tempmodel.reactions]
    mnxidlist = [i.notes['metanetx_id'] for i in tempmodel.reactions]
    for model in modellist:
        for rea in model.reactions:
            for rea2 in tempmodel.reactions:
                if rea.notes['metanetx_id'] == rea2.notes['metanetx_id']:
                    #index_i = mnxidlist.index(rea.notes['metanetx_id'])
                    #rea2 = tempmodel.reactions.get_by_id(ridlist[index_i])
                    if model.id not in rea2.notes['from']:
                        rea2.notes['from'].append(model.id)
                        #print(model.id)
                    if model.id not in ['Lreu_se_dna_cobra','Lreu_se_rast_cobra']:
                        if (rea2.gene_reaction_rule in rea.gene_reaction_rule) or (len(rea2.gene_reaction_rule)<len(rea.gene_reaction_rule)):
                            if 'Unknown' not in rea.gene_reaction_rule:
                                rea2.gene_reaction_rule = rea.gene_reaction_rule
print(len(tempmodel.reactions),len(tempmodel.genes))

#save gaplist!!!

cobra.io.save_json_model(tempmodel,'add_two.json')












