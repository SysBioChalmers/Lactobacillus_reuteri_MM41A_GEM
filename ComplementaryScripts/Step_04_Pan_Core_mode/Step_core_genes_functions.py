#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 6/2/20

"""Step_core_genes_functions.py
:description : script
:param : 
:returns: 
:rtype: 
"""
import os

import cobra


rea_list = ['ACHBS',
            'ACLS',
            # 'MALLAC',
            'NTPP10',
            'NTPP11',
            'NTPP9',
            'ORNt2',
            'PRAGSr',
            'SUCBZL',
            'TREt']


os.chdir('../../ComplementaryData/Step_03_Compare_Refine/')
print('----- loading model -----')
iHL622 = cobra.io.load_json_model('../../ModelFiles/iHL622.json')
iML1515 = cobra.io.read_sbml_model('../Initial_data/template_models/iML1515.xml')

for rea_i in rea_list:
    model = iML1515.copy()
    rea = model.reactions.get_by_id(rea_i)
    print(rea.name)

    # print(model.optimize())