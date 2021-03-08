#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 2019-09-04

"""Step_37_strain_simulate.py
:description : script
:param : 
:returns: 
:rtype: 
"""

import os

import cobra
import matplotlib
import pandas as pd

matplotlib.rc('font', family="Arial")

matplotlib.rcParams["font.family"] = 'Arial'  # 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Arial']
os.chdir('../../ComplementaryData/Step_04_Pan_Core_model/')
print('----- loading data -----')

name_list = ['100_23', '20_2', '3c6', 'ATCC53608', 'CF48_3A',
             'DSM200016', 'I5007', 'IRT', 'JCM1112',
             'LTH2584', 'LTH5448', 'MM2_3', 'MM4_1A', 'SD2112',
             'TD1', 'TMW1_112', 'TMW1_656', 'lpuph', 'mlc3', 'LR1',
             'LR10', 'LR11', 'LR12', 'LR13', 'LR14',
             'LR17', 'LR18', 'LR19', 'LR2', 'LR3',
             'LR4', 'LR6', 'LR7', 'LR8', 'LR9'
             ]
output_table = 'models_df_products_2.tsv'

models_df = pd.DataFrame(columns=['model_id', 'growth', 'reaset', 'metset', 'genset',
                                  'lac__L_c', 'ac_c', 'etoh_c', 'hista_c', 'fol_c', 'adeadocblp_c', 'ppoh_c', '13ppd_c',
                                  '3hppnl_c', 'dhap_c', 'mthgxl_c', '12ppd__R_c'])
for sp_name in name_list:

    # Lreu_tp_model = cobra.io.load_json_model(sp_name + 'Lreu_draft_3_refined_19-08-21.json')
    # Lreu_tp_model.solver = 'glpk'

    # keep_move_dic = {
    #     "MCOATA": "MACPMT",
    #     "DALTAL_LRE": "DALTAL",
    #     "KAS15": "kaasIII",
    #     "ALATRS": "ALATRS_1",
    #     "ASPTRS": "ASPTRS_1",
    #     "CYSTRS": "CYSTRS_1",
    #     "HISTRS": "HISTRS_1",
    #     "LEUTRS": "LEUTRS_1",
    #     "METTRS": "METTRS_1",
    #     "PHETRS": "PHETRS_1",
    #     "ARGTRS": "ARGTRS_1",
    #     "SERTRS": "SERTRS_1",
    #     "GLYTRS": "GLYTRS_1",
    #     "ILETRS": "ILETRS_1",
    #     "LYSTRS": "LYSTRS_1",
    #     "THRTRS": "THRTRS_1",
    #     "TRPTRS": "TRPTRS_1",
    #     "TYRTRS": "TYRTRS_1",
    #     "VALTRS": "VALTRS_1",
    #     "ASNTRS": "ASNTRS_1",
    #
    #     "EAR40x": "BTMAT1",
    #     "PRFGS": "PRFGS_1",
    #     "PPA2": "EPPP",
    #
    #     "LEUt2r": "LEUt2",
    #
    #     "ARBt2": "ARAB_Lt",
    #     "ARBt2r": "ARBt2",
    #     "SERt3": "SERt2r",
    #     "SERt2r": "SERt3",
    #     # "THRt2r":"THRt3",
    #
    #     "AGAT_LRE": "AGAT_LPL",
    #     "CLPNS_LRE": "CLPNS_LPL",
    #     "CPSS_LRE": "CPSS_LPL",
    #     # "CPSS_LRE":"CPSS_LPL2",
    #     "DAGGT_LRE": "DAGGT_LPL",
    #     "DALTAL_LPL": "DALTAL",
    #     "GLTAL_LPL": "GLTAL",
    #     "DASYN_LRE": "LPGS_LPL",
    #     "PGSA_LRE": "PGSA_LPL"
    # }
    #
    # for k, v in keep_move_dic.items():
    #     Lreu_tp_model = My_def.merge_model.merge_reactionsid(Lreu_tp_model, k, v)
    # cobra.io.save_json_model(Lreu_tp_model,sp_name + 'Lreu_draft_3_refined_v2.json')
    Lreu_tp_model = cobra.io.load_json_model(sp_name + 'Lreu_draft_3_refined_v2.json')

    # <get model information: table>

    reaset = set([i.id for i in Lreu_tp_model.reactions])
    growth_rate = Lreu_tp_model.optimize().objective_value

    # prod_list = ['lac__L_c','ac_c','etoh_c','hista_c','fol_c','adeadocbl_c','ppoh_c','13ppd_c']
    # prod_values= []

    prod_dict = {'lac__L_c': -1, 'ac_c': -1, 'etoh_c': -1, 'hista_c': -1,
                 'fol_c': -1, 'adeadocbl_c': -1, 'ppoh_c': -1, '13ppd_c': -1,
                 '3hppnl_c': -1, 'dhap_c': -1, 'mthgxl_c': -1, '12ppd__R_c': -1}
    # prod_dict= {'12ppd__R_c':-1,'mthgxl_c':-1,'dhap_c':-1,'hista_c':-1}
    prod_dict.update({'model_id': sp_name,
                      'reaset': reaset,
                      'growth': growth_rate})

    model = Lreu_tp_model.copy()
    model.reactions.get_by_id('EX_glc__D_e').bounds = (-20, 1000)
    model.reactions.get_by_id('EX_glyc_e').bounds = (-10, 1000)

    for prod in prod_dict.keys():
        model_ = model.copy()
        if 'hista_c' == prod:
            try:
                model_.reactions.get_by_id('HISDC')
                prod_dict[prod] = 1
            except:
                prod_dict[prod] = 0
        else:

            try:
                try:
                    model_.metabolites.get_by_id(prod)
                except:

                    continue
                rea = cobra.Reaction('Obj')
                model_.add_reaction(rea)
                model_.reactions.get_by_id('Obj').reaction = prod + ' --> '
                model_.objective = 'Obj'
                solution = model_.optimize()
                if solution.objective_value > 0.001:
                    prod_dict[prod] = solution.objective_value
                else:
                    prod_dict[prod] = 0
            except:
                prod_dict[prod] = 0

    models_df = models_df.append(prod_dict,
                                 ignore_index=True)

    models_df.to_csv(output_table, sep='\t')

# for i in models_df_.index:
#     sp_name = models_df_.iloc[i]['model_id']
#
#     locals()[sp_name + '_reaset'] = eval(models_df_.iloc[i]['reaset'])
#     locals()[sp_name + '_metset'] = eval(models_df_.iloc[i]['metset'])
#     locals()[sp_name + '_genset'] = eval(models_df_.iloc[i]['genset'])

# %%

her_name = [
    'LR1',
    'LR10', 'LR11', 'LR12', 'LR13', 'LR14',
    'LR17', 'LR18', 'LR19', 'LR2', 'LR3',
    'LR4', 'LR6', 'LR7', 'LR8', 'LR9'
]

omn_name = ['JCM1112', 'MM2_3', 'MM4_1A', 'CF48_3A', 'SD2112', 'I5007', 'ATCC53608', 'DSM200016', 'IRT', 'TD1', 'mlc3',
            '100_23', '20_2', '3c6', 'lpuph', ]
sour_name = ['LTH5448', 'TMW1_112', 'TMW1_656', 'LTH2584']

Lreu_d = {'sp_name': name_list, }
Lreu_df = pd.DataFrame(data=Lreu_d)
Lreu_df['group'] = 'her'

Lreu_df.loc[Lreu_df['sp_name'].isin(omn_name), ['group']] = 'omn'
Lreu_df.loc[Lreu_df['sp_name'].isin(sour_name), ['group']] = 'sou'
Lreu_df = Lreu_df.sort_values(by=['group', ])
models_df.to_csv(output_table, sep='\t')
# models_df_ = pd.read_csv('models_df_products.tsv', sep='\t', index_col=0)
