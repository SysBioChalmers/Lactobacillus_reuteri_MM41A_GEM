#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 2019-06-10

"""step_04_refine_model.py
:description : script
:param : 
:returns: 
:rtype: 
"""

import cobra
import os


os.chdir('../../ComplementaryData/Step_03_Compare_Refine/')

#%% load data
iNF517 = cobra.io.load_json_model('iNF517.json')
iML1515 = cobra.io.load_json_model('iML1515.json')
Lreu_merged = cobra.io.load_json_model('../Step_03_Compare_Refine/Lreu_merged_gapfiled.json')

Lreu_merged_initial = cobra.io.load_json_model('Lreu_merged.json')

ngset = set()
migset = set()
for i in Lreu_merged.genes:
    if 'MBLCL' in i.id:
        ngset.add(i.id)
    elif 'missing' in i.id:
        migset.add(i.id)

print('ng:' ,len(ngset))
print('migset:' ,len(migset))


infset = set()
ibtset = set()
imlset = set()
for i in Lreu_merged_initial.reactions:
    if 'iNF517' in i.notes['from'] or 'Lreu_from_iNF517' in i.notes['from']:
        infset.add(i.id)
    if 'iBT721' in i.notes['from'] or 'Lreu_from_iBT721' in i.notes['from']:
        ibtset.add(i.id)
    if 'iML1515' in i.notes['from'] or 'Lreu_from_iML1515' in i.notes['from']:
        imlset.add(i.id)
print('517',len(infset))
print('721',len(ibtset - infset))
print('1515',len(imlset - infset - ibtset))




exset = set()
transet = set()
biomassset = set()

for rea in Lreu_merged.reactions:
    if ('EX_' in rea.id):
        #exchange reactions
        exset.add(rea.id)
    elif ('_c' in rea.reaction) and ('_e' in rea.reaction):
        transet.add(rea.id)

    elif ('_LRE' in rea.id) or ('LRE_c' in rea.reaction):
        biomassset.add(rea.id)

print('exset',len(exset))
print('transet',len(transet))
print('init',len(Lreu_merged.reactions) - len(transet| transet)-len(exset))
print('biomassset',len(biomassset))



exmset = set()
cmnset = set()

for rea in Lreu_merged.metabolites:
    if rea.id.endswith('_e'):
        #exchange reactions
        exmset.add(rea.id)
    elif rea.id.endswith('_c'):
        cmnset.add(rea.id)

print('exmset',len(exmset))
print('cmnset',len(cmnset))

# %%
for rea in Lreu_merged.reactions:
    if 'EX' in rea.id:
        if rea.lower_bound<=0 and rea.upper_bound <= 0:
            rea.upper_bound = 0.0
        elif rea.lower_bound >= 0 and rea.upper_bound >= 0:
            rea.lower_bound = 0.0

for rea in iNF517.reactions:
    if 'EX' in rea.id:
        if rea.lower_bound<=0 and rea.upper_bound <= 0:
            rea.upper_bound = 0.0
        elif rea.lower_bound >= 0 and rea.upper_bound >= 0:
            rea.lower_bound = 0.0



Lreu_merged.objective = "BIOMASS_LRE"
print('Biomass:',Lreu_merged.optimize())

iNF517.objective = "BIOMASS_LRE"
print('Biomass:',iNF517.optimize())

# %% simulation
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
glyc_reaid = 'EX_glyc_e'
glc_reaid = 'EX_glc__D_e'

aadic = {}
for i in aarealist:
    rea1 = Lreu_merged.reactions.get_by_id(i)
    rea2 = iNF517.reactions.get_by_id(i)
    bounds1 = rea1.bounds
    bounds2 = rea2.bounds

    rea1.bounds = (0.0,10)
    rea2.bounds = (0.0,10)
    #print('----- %s -----'%i )
    aadic[i] = [Lreu_merged.optimize().objective_value, iNF517.optimize().objective_value]

    rea1.bounds = bounds1
    rea2.bounds = bounds2

print (Lreu_merged.optimize().objective_value, iNF517.optimize().objective_value)
b1 = Lreu_merged.reactions.get_by_id('EX_glyc_e').lower_bound
b2 = Lreu_merged.reactions.get_by_id('EX_glc__D_e').lower_bound

Lreu_merged.reactions.get_by_id('EX_glyc_e').lower_bound = - 2.2
Lreu_merged.reactions.get_by_id('EX_glc__D_e').lower_bound = - 1.52

iNF517.reactions.get_by_id('EX_glyc_e').lower_bound = - 2.2
iNF517.reactions.get_by_id('EX_glc__D_e').lower_bound = - 1.52


print (Lreu_merged.optimize().objective_value, iNF517.optimize().objective_value)


Lreu_merged.reactions.get_by_id('EX_glyc_e').lower_bound = b1
Lreu_merged.reactions.get_by_id('EX_glc__D_e').lower_bound = b2

iNF517.reactions.get_by_id('EX_glyc_e').lower_bound = b1
iNF517.reactions.get_by_id('EX_glc__D_e').lower_bound = b2




























