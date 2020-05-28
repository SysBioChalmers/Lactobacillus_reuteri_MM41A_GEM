#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 4/22/19

"""Step_refine_pipeline_part03_amino_acids.py
:description : script to refine amino acids
:param : 
:returns: 
:rtype: 
"""

import os

import cobra

import My_def

os.chdir('../../ComplementaryData/Step_03_Compare_Refine/')

Lreu_draft_3_refined = cobra.io.load_json_model('Lreu_draft_3_refined_part02.json')
Lreuteri_530 = cobra.io.load_json_model('../Step_02_DraftModels/Template/template_models/Lreuteri_530_standlized.json')
iNF517 = cobra.io.load_json_model('../Step_02_DraftModels/Template/template_models/iNF517_standlized.json')
iML1515 = cobra.io.load_json_model('../Step_02_DraftModels/Template/template_models/iML1515_standlized.json')

# %% <check growth>

print(Lreu_draft_3_refined.optimize().objective_value)
print(Lreuteri_530.optimize().objective_value)
print(iNF517.optimize().objective_value)
print(iML1515.optimize().objective_value)
# Lreu_draft_3_refined.reactions.get_by_id('PFK').bounds = (0,1000)


# %% <check mediume and basic carbon source and ATP>

print(Lreu_draft_3_refined.medium)
print(Lreuteri_530.medium)
print(iNF517.medium)
print(iML1515.medium)

Lreu_draft_3_refined.reactions.get_by_id('EX_nh4_e').bounds = (-1000, 1000)
Lreu_draft_3_refined.reactions.get_by_id('PFK').bounds = (0, 1000)  # NOTE： Lreuteri_530 is（-1000，2）

# Lreu_draft_3_refined.add_reaction(Lreuteri_530.reactions.get_by_id('GLCt'))
# Lreu_draft_3_refined.reactions.get_by_id('EX_nh3_e').bounds = (-1000,1000)
# Lreu_draft_3_refined.add_reaction(Lreuteri_530.reactions.get_by_id('NH3c'))


# other carbon:
# Lreu_draft_3_refined.reactions.get_by_id('EX_cys__L_e').bounds = (0,1000)
# Lreu_draft_3_refined.reactions.get_by_id('EX_gly_e').bounds = (0,1000)
# Lreu_draft_3_refined.reactions.get_by_id('EX_ala__L_e').bounds = (0,1000)
# Lreu_draft_3_refined.reactions.get_by_id('EX_leu__L_e').bounds = (0,1000)
# Lreu_draft_3_refined.reactions.get_by_id('EX_ile__L_e').bounds = (0,1000)
# Lreu_draft_3_refined.reactions.get_by_id('EX_thr__L_e').bounds = (0,1000)


# uselsee:
# Lreu_draft_3_refined.reactions.get_by_id('2H3MBt').remove_from_model()
# Lreu_draft_3_refined.reactions.get_by_id('2H3MPt').remove_from_model()

# Lreu_draft_3_refined.reactions.get_by_id('EX_2h3mp_e').bounds = (0,0)
# Lreu_draft_3_refined.reactions.get_by_id('EX_2h3mb_e').bounds = (0,0)
# Lreu_draft_3_refined.reactions.get_by_id('EX_2hxic__L_e').bounds = (0,0)
# Lreu_draft_3_refined.reactions.get_by_id('EX_ade_e').bounds = (0,0)
# Lreu_draft_3_refined.reactions.get_by_id('EX_o2_e').bounds = (-1000,1000)


for model_i in [Lreu_draft_3_refined, Lreuteri_530, iML1515]:  # iNF517
    model = model_i.copy()

    try:
        # other carbon:
        model.reactions.get_by_id('EX_cys__L_e').bounds = (0, 1000)
        model.reactions.get_by_id('EX_gly_e').bounds = (0, 1000)
        model.reactions.get_by_id('EX_ala__L_e').bounds = (0, 1000)
        model.reactions.get_by_id('EX_leu__L_e').bounds = (0, 1000)
        model.reactions.get_by_id('EX_ile__L_e').bounds = (0, 1000)
        model.reactions.get_by_id('EX_thr__L_e').bounds = (0, 1000)

        model.reactions.get_by_id('EX_arg__L_e').bounds = (0, 1000)
        model.reactions.get_by_id('EX_asn__L_e').bounds = (0, 1000)
        model.reactions.get_by_id('EX_asp__L_e').bounds = (0, 1000)
        model.reactions.get_by_id('EX_glu__L_e').bounds = (0, 1000)
        model.reactions.get_by_id('EX_glu__L_e').bounds = (0, 1000)
        model.reactions.get_by_id('EX_glu__L_e').bounds = (0, 1000)
        model.reactions.get_by_id('EX_ser__L_e').bounds = (0, 1000)

        # other limations:
        model.reactions.get_by_id('EX_etoh_e').bounds = (0, 1000)
        model.reactions.get_by_id('EX_lac__L_e').bounds = (0, 1000)
        model.reactions.get_by_id('EX_lac__L_e').bounds = (0, 1000)
        model.reactions.get_by_id('EX_o2_e').bounds = (0, 1000)

    except:
        pass

    model.reactions.get_by_id('ATPM').bounds = (0, 1000)
    model.reactions.get_by_id('EX_glc__D_e').bounds = (-10, 1000)

    model.objective = 'EX_co2_e'

    # model.objective = 'ATPM'
    # rea = cobra.Reaction('NADHM')
    # model.add_reaction(rea)
    # model.reactions.get_by_id('NADHM').reaction = ' nadh_c--> nad_c'
    # model.objective = 'NADHM'

    solution = model.optimize()
    print(solution)

    solution = cobra.flux_analysis.pfba(model)  # model.optimize()
    My_def.io_file.solution2txt(solution, model, model.id + '_temp_flux.txt')

# %% amino acids

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

unessential = [
    'EX_ala__L_e',
    'EX_asp__L_e',
    'EX_cys__L_e',
    'EX_gly_e',
    'EX_ile__L_e',
    'EX_lys__L_e',
    'EX_pro__L_e',
    'EX_ser__L_e', ]

essential = [
    'EX_arg__L_e',
    'EX_asn__L_e',
    'EX_gln__L_e',
    'EX_glu__L_e',
    'EX_his__L_e',
    'EX_leu__L_e',
    'EX_met__L_e',
    'EX_phe__L_e',
    'EX_thr__L_e',
    'EX_trp__L_e',
    'EX_tyr__L_e',
    'EX_val__L_e']

# EX_glyc_e: Glycerol

Lreu_draft_3_refined.objective = 'BIOMASS'
print('initial biomass: \t', Lreu_draft_3_refined.optimize())

Lreu_draft_3_refined.reactions.get_by_id('SERt2r').reaction = Lreuteri_530.reactions.get_by_id('SERt2r').reaction
Lreu_draft_3_refined.reactions.get_by_id('araphe1').bounds = (0, 0)
Lreu_draft_3_refined.reactions.get_by_id('aratyr1').bounds = (0.0, 0.0)

Lreu_draft_3_refined.reactions.get_by_id('EX_ala__L_e').bounds = (-2.69, 1000)
Lreu_draft_3_refined.reactions.get_by_id('EX_asp__L_e').bounds = (-3.16, 1000)
Lreu_draft_3_refined.reactions.get_by_id('EX_cys__L_e').bounds = (-0.83, 1000)
Lreu_draft_3_refined.reactions.get_by_id('EX_gly_e').bounds = (-2.33, 1000)
Lreu_draft_3_refined.reactions.get_by_id('EX_ile__L_e').bounds = (-1.60, 1000)
Lreu_draft_3_refined.reactions.get_by_id('EX_lys__L_e').bounds = (-2.68, 1000)
Lreu_draft_3_refined.reactions.get_by_id('EX_pro__L_e').bounds = (-5.86, 1000)
Lreu_draft_3_refined.reactions.get_by_id('EX_ser__L_e').bounds = (-3.24, 1000)

Lreu_draft_3_refined.reactions.get_by_id('EX_gln__L_e').bounds = (-1, 1000)
Lreu_draft_3_refined.reactions.get_by_id('EX_thr__L_e').bounds = (-1, 1000)
Lreu_draft_3_refined.reactions.get_by_id('EX_leu__L_e').bounds = (-1, 1000)

Lreu_draft_3_refined.objective = 'BIOMASS'
print('ref amino acid biomass: \t', Lreu_draft_3_refined.optimize())

# ile
Lreu_draft_3_refined.add_reaction(iNF517.reactions.get_by_id('DHAD2'))
Lreu_draft_3_refined.add_reaction(iNF517.reactions.get_by_id('KARA2'))

# Cys:
Lreu_draft_3_refined.add_reaction(iML1515.reactions.get_by_id('SERAT'))
# Lreu_draft_3_refined.add_reaction(iML1515.reactions.get_by_id('CYSS'))
# Lreu_draft_3_refined.add_reaction(iML1515.reactions.get_by_id('SADT2'))
# Lreu_draft_3_refined.add_reaction(iML1515.reactions.get_by_id('ADSK'))
# Lreu_draft_3_refined.add_reaction(iML1515.reactions.get_by_id('PAPSR'))
# Lreu_draft_3_refined.add_reaction(iML1515.reactions.get_by_id('BPNT'))

print('refined model biomass: \t', Lreu_draft_3_refined.optimize())

# thr based on iML1515
Lreu_draft_3_refined.reactions.get_by_id('THRA').bounds = (0, 1000)

# thr and tyr based on Lreuteri_530
Lreu_draft_3_refined.reactions.get_by_id('aratry1').bounds = (0, 0)
Lreu_draft_3_refined.reactions.get_by_id('araphe3').bounds = (0, 0)
# pro
# Lreu_draft_3_refined.add_reaction(iML1515.reactions.get_by_id('CYSS'))
# Lreu_draft_3_refined.add_reaction(iML1515.reactions.get_by_id('SADT2'))


# ? fale gap
# Lreu_draft_3_refined.reactions.get_by_id('TRPS1').bounds = (0, 0)
# Lreu_draft_3_refined.reactions.get_by_id('PPND').bounds = (0, 0)
# Lreu_draft_3_refined.reactions.get_by_id('TRPS1').remove_from_model()
# Lreu_draft_3_refined.reactions.get_by_id('PPND').remove_from_model()

# %%
# Lreu_draft_3_refined.reactions.get_by_id('ASPt2r').bounds = (-1000, 1000)
# Lreu_draft_3_refined.reactions.get_by_id('CYSt2r').bounds = (-1000, 1000)
# Lreu_draft_3_refined.reactions.get_by_id('SERt2r').bounds = (-1000, 1000)
# Lreu_draft_3_refined.reactions.get_by_id('ASNt2r').bounds = (-1000, 1000)
# Lreu_draft_3_refined.reactions.get_by_id('ARGt2r').bounds = (-1000, 1000)
# Lreu_draft_3_refined.reactions.get_by_id('GLNabc').bounds = (-1000, 1000)


aarealist = essential
aarealist = unessential
aarealist = unessential + essential

# Lreu_draft_3_refined.reactions.get_by_id('EX_gln__L_e').bounds = (-1, 1000)
# Lreu_draft_3_refined.reactions.get_by_id('EX_glu__L_e').bounds = (-1, 1000)
#
# Lreu_draft_3_refined.reactions.get_by_id('EX_asp__L_e').bounds = (0, 1000)
# Lreu_draft_3_refined.reactions.get_by_id('EX_arg__L_e').bounds = (-1, 1000)
# Lreu_draft_3_refined.reactions.get_by_id('EX_asn__L_e').bounds = (-1, 1000)

# Lreu_draft_3_refined.reactions.get_by_id('G5SADs').bounds = (0, 0)
# Lreu_draft_3_refined.reactions.get_by_id('HSDy').bounds = (0, 0)

aadic = {}
for aa in aarealist:
    aadic[aa] = []
    for model_i in [Lreu_draft_3_refined, Lreuteri_530, ]:  # ,iNF517,iML1515
        model = model_i.copy()
        model.objective = 'BIOMASS'  # aa / 'BIOMASS'
        rea = model.reactions.get_by_id(aa)
        rea.bounds = (0, 1000)
        solution = model.optimize()  # cobra.flux_analysis.pfba(model)
        aadic[aa].append(solution)
        solution = cobra.flux_analysis.pfba(model)  # cobra.flux_analysis.pfba(model)
        # My_def.io_file.solution2txt(solution, model, model.id + aa + '_temp_flux.txt')
    print(aa, aadic[aa])

Lreu_draft_3_refined.objective = 'BIOMASS'
print('ref amino acid biomass: \t', Lreu_draft_3_refined.optimize().objective_value)

cobra.io.save_json_model(Lreu_draft_3_refined, 'Lreu_draft_3_refined_part03.json')
My_def.io_file.model2txt(Lreu_draft_3_refined, 'Lreu_draft_3_refined_part03.txt', sort=True)

cobra.io.write_sbml_model(Lreu_draft_3_refined, 'Lreu_draft_3_refined_part03.xml')

comd = ' memote report snapshot --filename "Lreu_draft_3_refined_part03.html" Lreu_draft_3_refined_part03.xml'
os.system(comd)