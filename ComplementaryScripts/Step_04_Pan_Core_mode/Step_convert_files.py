#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 3/8/21

"""Step_convert_files.py
:description : script
:param : 
:returns: 
:rtype: 
"""

import os

import cobra

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
for sp_name in name_list:
    Lreu_tp_model = cobra.io.load_json_model(sp_name + 'Lreu_draft_3_refined_v2.json')

    # cobra.io.write_sbml_model(Lreu_tp_model, '../../ModelFiles/pan_core_GEMs/Lreu_' + sp_name + '.xml')
    cobra.io.save_json_model(Lreu_tp_model, '../../ModelFiles/pan_core_GEMs/Lreu_' + sp_name + '.json', sort='True')
