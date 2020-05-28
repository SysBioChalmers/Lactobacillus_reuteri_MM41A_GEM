#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 4/22/20

"""Step_refine_pipeline_part03_amino_acids.py
:description : script to refine important metabolites
:param : metabolites list: [lac acet etoh 1,3-propanediol,reuterin,hista_HISDC]
:returns: 
:rtype:

"""

import os

import cobra

import My_def

os.chdir('../../ComplementaryData/Step_03_Compare_Refine/')

# %%
Lreu_draft_3_refined = cobra.io.load_json_model('Lreu_draft_3_refined_part01.json')
Lreu_draft_3_refined.id = 'Lreu_draft_3_refined'
Lreuteri_530 = cobra.io.load_json_model('../Step_02_DraftModels/Template/template_models/Lreuteri_530_standlized.json')
iNF517 = cobra.io.load_json_model('../Step_02_DraftModels/Template/template_models/iNF517_standlized.json')
iML1515 = cobra.io.load_json_model('../Step_02_DraftModels/Template/template_models/iML1515_standlized.json')
iML1515.reactions.get_by_id('EX_glc__D_e').bounds = (-25, 1000)

# %% <lac: >
met_i = Lreu_draft_3_refined.metabolites.get_by_id('lac__L_e')
print('\n', met_i.id, met_i.name, met_i.formula, met_i.annotation)
Lreu_draft_3_refined.reactions.get_by_id('EX_lac__L_e').bounds = (0, 1000)

objective_rea = 'EX_lac__L_e'
Lreu_draft_3_refined.objective = objective_rea
print('Lreu opt value:  ', Lreu_draft_3_refined.optimize().objective_value)
# experiment data: Y LACt/GLU = 0.469 0.682 0.076

for model_i in [Lreu_draft_3_refined, Lreuteri_530, iML1515, ]:  # iNF517
    model = model_i.copy()

    model.reactions.get_by_id('EX_glc__D_e').bounds = (-25, 1000)
    model.reactions.get_by_id(objective_rea).bounds = (0, 1000)
    model.reactions.get_by_id('ATPM').bounds = (0, 1000)
    model.reactions.get_by_id('PFK').bounds = (0, 1000)  # NOTE： Lreuteri_530 is（-1000，2）

    try:
        # close other carbon:
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

    model.objective = objective_rea
    solution = model.optimize()
    print(model.id, solution.objective_value)
    obj_value = solution.objective_value

    Y_g_g = obj_value * model.metabolites.get_by_id('lac__L_c').formula_weight \
            / (25 * model.metabolites.get_by_id('glc__D_e').formula_weight)

    Y_c = obj_value * 3 / (25 * 6)

    print('Y_g_g: ', Y_g_g, 'T_c: ', Y_c)

    solution = cobra.flux_analysis.pfba(model)  # model.optimize()
    My_def.io_file.solution2txt(solution, model, model.id + '_temp_flux.txt')

# %% <acet:> TODO:
met_i = Lreu_draft_3_refined.metabolites.get_by_id('ac_e')
print('\n', met_i.id, met_i.name, met_i.formula, met_i.annotation)
Lreu_draft_3_refined.reactions.get_by_id('EX_ac_e').bounds = (0, 1000)
Lreu_draft_3_refined.add_reaction(iML1515.reactions.get_by_id('SADT2'))
Lreu_draft_3_refined.add_reaction(iML1515.reactions.get_by_id('ADSK'))
Lreu_draft_3_refined.add_reaction(iML1515.reactions.get_by_id('PAPSR'))
Lreu_draft_3_refined.add_reaction(iML1515.reactions.get_by_id('BPNT'))

objective_rea = 'EX_ac_e'  # 'EX_co2_e'#'EX_for_e'##,
Lreu_draft_3_refined.objective = objective_rea
print('Lreu opt value:  ', Lreu_draft_3_refined.optimize().objective_value)
# experiment data: Y ac/GLU = 0.2170-0.33

for model_i in [Lreu_draft_3_refined, Lreuteri_530, iML1515, ]:  # iNF517
    model = model_i.copy()

    model.reactions.get_by_id('EX_glc__D_e').bounds = (-25, 1000)
    model.reactions.get_by_id(objective_rea).bounds = (0, 1000)
    model.reactions.get_by_id('ATPM').bounds = (0, 1000)
    model.reactions.get_by_id('PFK').bounds = (0, 1000)  # NOTE： Lreuteri_530 is（-1000，2）

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
        model.reactions.get_by_id('EX_for_e').bounds = (0, 1000)
        model.reactions.get_by_id('EX_etoh_e').bounds = (0, 1000)
        model.reactions.get_by_id('EX_lac__L_e').bounds = (0, 1000)
        model.reactions.get_by_id('EX_lac__L_e').bounds = (0, 1000)
        model.reactions.get_by_id('EX_o2_e').bounds = (0, 1000)

    except:
        pass
    # rea = cobra.Reaction('NADHM')
    # model.add_reaction(rea)
    # model.reactions.get_by_id('NADHM').reaction = 'nadh_c + h_c --> nad_c'
    # model.objective = 'NADHM'

    model.objective = objective_rea
    solution = model.optimize()
    print(model.id, solution.objective_value)
    obj_value = solution.objective_value

    Y_g_g = obj_value * model.metabolites.get_by_id('ac_e').formula_weight \
            / (25 * model.metabolites.get_by_id('glc__D_e').formula_weight)

    Y_c = obj_value * 2 / (25 * 6)

    print('Y_g_g: ', Y_g_g, 'T_c: ', Y_c)

    solution = cobra.flux_analysis.pfba(model)  # model.optimize()
    My_def.io_file.solution2txt(solution, model, model.id + '_temp_flux.txt')

# %% <etoh:>
met_i = Lreu_draft_3_refined.metabolites.get_by_id('etoh_e')
print('\n', met_i.id, met_i.name, met_i.formula, met_i.annotation)
Lreu_draft_3_refined.reactions.get_by_id('EX_etoh_e').bounds = (0, 1000)

objective_rea = 'EX_etoh_e'
Lreu_draft_3_refined.objective = objective_rea
print('Lreu opt value:  ', Lreu_draft_3_refined.optimize().objective_value)
# experiment data: Y ac/GLU = 0.2170-0.33

for model_i in [Lreu_draft_3_refined, Lreuteri_530, iML1515, ]:  # iNF517
    model = model_i.copy()

    model.reactions.get_by_id('EX_glc__D_e').bounds = (-25, 1000)
    model.reactions.get_by_id(objective_rea).bounds = (0, 1000)
    model.reactions.get_by_id('ATPM').bounds = (0, 1000)
    model.reactions.get_by_id('PFK').bounds = (0, 1000)  # NOTE： Lreuteri_530 is（-1000，2）

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
    # rea = cobra.Reaction('NADHM')
    # model.add_reaction(rea)
    # model.reactions.get_by_id('NADHM').reaction = 'nadh_c + h_c --> nad_c'
    # model.objective = 'NADHM'

    model.objective = objective_rea
    solution = model.optimize()
    print(model.id, solution.objective_value)
    obj_value = solution.objective_value

    Y_g_g = obj_value * model.metabolites.get_by_id('etoh_e').formula_weight \
            / (25 * model.metabolites.get_by_id('glc__D_e').formula_weight)

    Y_c = obj_value * 2 / (25 * 6)

    print('Y_g_g: ', Y_g_g, 'T_c: ', Y_c)

    solution = cobra.flux_analysis.pfba(model)  # model.optimize()
    My_def.io_file.solution2txt(solution, model, model.id + '_temp_flux.txt')

# %%<1-propanol: ppoh> TODO:
met_i = Lreu_draft_3_refined.metabolites.get_by_id('ppoh_c')
print('\n', met_i.id, met_i.name, met_i.formula, met_i.annotation)

Lreu_draft_3_refined.add_reaction(iML1515.reactions.get_by_id('MGSA'))
Lreu_draft_3_refined.add_reaction(Lreuteri_530.reactions.get_by_id('PPOHt'))
Lreu_draft_3_refined.add_reaction(Lreuteri_530.reactions.get_by_id('EX_ppoh_e'))
Lreu_draft_3_refined.reactions.get_by_id('EX_ppoh_e').bounds = (0, 1000)

objective_rea = 'EX_ppoh_e'
Lreu_draft_3_refined.objective = objective_rea
print('Lreu opt value:  ', Lreu_draft_3_refined.optimize().objective_value)
# experiment data: Y ac/GLU = 0.2170-0.33

for model_i in [Lreu_draft_3_refined, Lreuteri_530, ]:  # iNF517iML1515,
    model = model_i.copy()

    model.reactions.get_by_id('EX_glc__D_e').bounds = (-25, 1000)
    # model.reactions.get_by_id('EX_glyc_e').bounds = (-10, 1000)
    model.reactions.get_by_id(objective_rea).bounds = (0, 1000)
    model.reactions.get_by_id('ATPM').bounds = (0, 1000)
    model.reactions.get_by_id('PFK').bounds = (0, 1000)  # NOTE： Lreuteri_530 is（-1000，2）

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
    # rea = cobra.Reaction('NADHM')
    # model.add_reaction(rea)
    # model.reactions.get_by_id('NADHM').reaction = 'nadh_c + h_c --> nad_c'
    # model.objective = 'NADHM'

    model.objective = objective_rea
    solution = model.optimize()
    print(model.id, solution.objective_value)
    obj_value = solution.objective_value

    Y_g_g = obj_value * model.metabolites.get_by_id('ppoh_c').formula_weight \
            / (25 * model.metabolites.get_by_id('glc__D_e').formula_weight)

    Y_c = obj_value * 3 / (25 * 6)

    print('Y_g_g: ', Y_g_g, 'T_c: ', Y_c)

    solution = cobra.flux_analysis.pfba(model)  # model.optimize()
    My_def.io_file.solution2txt(solution, model, model.id + '_temp_flux.txt')

# %% <1,3-propanediol ~> TODO:
met_i = Lreu_draft_3_refined.metabolites.get_by_id('13ppd_c')
print('\n', met_i.id, met_i.name, met_i.formula, met_i.annotation)

# Lreu_draft_3_refined.add_reaction(Lreuteri_530.reactions.get_by_id('PPN13D'))
Lreu_draft_3_refined.add_reaction(Lreuteri_530.reactions.get_by_id('PPDt1'))
Lreu_draft_3_refined.add_reaction(Lreuteri_530.reactions.get_by_id('EX_13ppd_e'))

# Lreu_draft_3_refined.reactions.get_by_id('EX_glyc_e').bounds = (0,1000)
Lreu_draft_3_refined.reactions.get_by_id('EX_13ppd_e').bounds = (0, 1000)

objective_rea = 'EX_13ppd_e'
Lreu_draft_3_refined.objective = objective_rea
print('Lreu opt value:  ', Lreu_draft_3_refined.optimize().objective_value)
# Y  = Lreu_draft_3_refined.optimize().objective_value * Lreu_draft_3_refined.metabolites.get_by_id('13ppd_e').formula_weight\
#      /(10*Lreu_draft_3_refined.metabolites.get_by_id('glc__D_e').formula_weight )
# solution = Lreu_draft_3_refined.optimize()
# print(Y)
# My_def.io_file.solution2txt(solution,Lreu_draft_3_refined,'Lreu_draft_3_refined_temp_flux.txt')

for model_i in [Lreu_draft_3_refined, Lreuteri_530, ]:  # iNF517,iML1515,
    model = model_i.copy()

    model.reactions.get_by_id('EX_glc__D_e').bounds = (-25, 1000)
    model.reactions.get_by_id(objective_rea).bounds = (0, 1000)
    model.reactions.get_by_id('ATPM').bounds = (0, 1000)
    model.reactions.get_by_id('PFK').bounds = (0, 1000)  # NOTE： Lreuteri_530 is（-1000，2）

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
    # rea = cobra.Reaction('NADHM')
    # model.add_reaction(rea)
    # model.reactions.get_by_id('NADHM').reaction = 'nadh_c + h_c --> nad_c'
    # model.objective = 'NADHM'

    model.objective = objective_rea
    solution = model.optimize()
    print(model.id, solution.objective_value)
    obj_value = solution.objective_value

    Y_g_g = obj_value * model.metabolites.get_by_id('13ppd_e').formula_weight \
            / (25 * model.metabolites.get_by_id('glc__D_e').formula_weight)

    Y_c = obj_value * 3 / (25 * 6)

    print('Y_g_g: ', Y_g_g, 'T_c: ', Y_c)

    solution = cobra.flux_analysis.pfba(model)  # model.optimize()
    My_def.io_file.solution2txt(solution, model, model.id + '_temp_flux.txt')

# %%<reuterin : 3 hydroxypropionaldehyde, 3-HPA,3 hydroxypropanal, 3hppnl> : TODO: 3hppnl
met_i = Lreu_draft_3_refined.metabolites.get_by_id('3hpp_e')
print('\n', met_i.id, met_i.name, met_i.formula, met_i.annotation)
Lreu_draft_3_refined.add_reaction(Lreuteri_530.reactions.get_by_id('3HPPt'))
Lreu_draft_3_refined.add_reaction(Lreuteri_530.reactions.get_by_id('EX_3hpp_e'))

objective_rea = 'EX_3hpp_e'
Lreu_draft_3_refined.objective = objective_rea
print('Lreu opt value:  ', Lreu_draft_3_refined.optimize().objective_value)

for model_i in [Lreu_draft_3_refined, Lreuteri_530, ]:  # iNF517,iML1515,
    model = model_i.copy()

    model.reactions.get_by_id('EX_glc__D_e').bounds = (-25, 1000)
    model.reactions.get_by_id(objective_rea).bounds = (0, 1000)
    model.reactions.get_by_id('ATPM').bounds = (0, 1000)
    model.reactions.get_by_id('PFK').bounds = (0, 1000)  # NOTE： Lreuteri_530 is（-1000，2）

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
    # rea = cobra.Reaction('NADHM')
    # model.add_reaction(rea)
    # model.reactions.get_by_id('NADHM').reaction = 'nadh_c + h_c --> nad_c'
    # model.objective = 'NADHM'

    model.objective = objective_rea
    # iML1515: 13PPDH2
    solution = model.optimize()
    print(model.id, solution.objective_value)
    obj_value = solution.objective_value

    Y_g_g = obj_value * model.metabolites.get_by_id('3hpp_e').formula_weight \
            / (25 * model.metabolites.get_by_id('glc__D_e').formula_weight)

    Y_c = obj_value * 3 / (25 * 6)

    print('Y_g_g: ', Y_g_g, 'T_c: ', Y_c)

    solution = cobra.flux_analysis.pfba(model)  # model.optimize()
    My_def.io_file.solution2txt(solution, model, model.id + '_temp_flux.txt')

# %% <Histamine hista_c>
met_i = Lreu_draft_3_refined.metabolites.get_by_id('hista_c')
print('\n', met_i.id, met_i.name, met_i.formula, met_i.annotation)

# Lreu_draft_3_refined.add_reaction(Lreuteri_530.reactions.get_by_id('HISDC'))
Lreu_draft_3_refined.add_reaction(Lreuteri_530.reactions.get_by_id('HISTAap'))
rea = cobra.Reaction('EX_hista_e')
Lreu_draft_3_refined.add_reaction(rea)
Lreu_draft_3_refined.reactions.get_by_id('EX_hista_e').reaction = 'hista_e --> '

objective_rea = 'EX_hista_e'
Lreu_draft_3_refined.objective = objective_rea
Lreu_draft_3_refined.reactions.get_by_id('EX_his__L_e').bounds = (-5, 1000)
print('Lreu opt  value:  ', Lreu_draft_3_refined.optimize().objective_value)

# %% <vitamin B12: cobalamin> TODO: Adenosylcobalamin, adocbl TODO :gap!!!!

# objectiverea = 'EX_adocbl_e'
# Lreu_draft_3_refined.objective = objectiverea
# print('Lreu opt  value:  ',Lreu_draft_3_refined.optimize().objective_value)
#
# Lreuteri_530.objective = objectiverea
# Lreuteri_530.optimize()

# NOTE in adocbl named adeadocbl, the pathway is not reliable.
# Lreu_draft_3_refined.add_reaction(iML1515.reactions.get_by_id('ADOCBLS'))
# Lreu_draft_3_refined.add_reaction(iML1515.reactions.get_by_id('RZ5PP'))
# Lreu_draft_3_refined.add_reaction(iML1515.reactions.get_by_id('NNDMBRT'))
#
# # Lreu_draft_3_refined.add_reaction(Lreuteri_530.reactions.get_by_id('HISTAap'))
#
# rea1 = cobra.Reaction('ADOCBLt')
# rea2 = cobra.Reaction('EX_adocbl_e')
# Lreu_draft_3_refined.add_reaction(rea1)
# Lreu_draft_3_refined.add_reaction(rea2)
# Lreu_draft_3_refined.reactions.get_by_id('ADOCBLt').reaction = 'adocbl_c --> adocbl_e'
# Lreu_draft_3_refined.reactions.get_by_id('EX_adocbl_e').reaction = 'adocbl_e --> '

# objectiverea = 'ADOCBLS'#'EX_adocbl_e'
# Lreu_draft_3_refined.objective = objectiverea
# print('Lreu opt  value:  ',Lreu_draft_3_refined.optimize().objective_value)
# iML1515.objective = 'BIOMASS_Ec_iML1515_WT_75p37M'
# iML1515.optimize()
# # TODO adeadocbl_c --> c


# %% <vitamin B9 Folate > TODO: fol_c
met_i = Lreu_draft_3_refined.metabolites.get_by_id('fol_e')
print('\n', met_i.id, met_i.name, met_i.formula, met_i.annotation)

objective_rea = 'EX_fol_e'  # 'EX_adocbl_e'
Lreu_draft_3_refined.objective = objective_rea
print('Lreu opt value:  ', Lreu_draft_3_refined.optimize().objective_value)

# %% <EPS: exopolysaccharide > NOTE: no
# %% <other SCFAs>
# %%<others:>
# <Mannitol mnl_e>
objective_rea = 'EX_mnl_e'  # 'EX_adocbl_e'
Lreu_draft_3_refined.objective = objective_rea
print('Lreu opt value:  ', Lreu_draft_3_refined.optimize().objective_value)

# %% <output files>
for i in Lreu_draft_3_refined.metabolites:
    if i.compartment not in ['c', 'e']:
        i.compartment = i.id.split('_')[-1]

Lreu_draft_3_refined.reactions.get_by_id('BIOMASS').reaction
Lreu_draft_3_refined.objective = 'BIOMASS'
print('Lreu opt biomass value:  ', Lreu_draft_3_refined.optimize().objective_value)

cobra.io.save_json_model(Lreu_draft_3_refined, 'Lreu_draft_3_refined_part02.json')
My_def.io_file.model2txt(Lreu_draft_3_refined, 'Lreu_draft_3_refined_part02.txt', sort=True)
cobra.io.write_sbml_model(Lreu_draft_3_refined, 'Lreu_draft_3_refined_part02.xml')

comd = ' memote report snapshot --filename "Lreu_draft_3_refined_part02.html" Lreu_draft_3_refined_part02.xml'
os.system(comd)
