#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 2019-08-21

"""Step_37strains_CarveMe.py
:description : script
:param : 
:returns: 
:rtype: 
"""


import os
import pandas as pd
import re
import My_def
import cobra
from importlib import  reload
reload(My_def)

pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
name_list = ['100_23','20_2','3c6','ATCC53608','CF48_3A',
             'DSM200016','I5007','IRT','JCM1112',
             'LTH2584','LTH5448','MM2_3','MM4_1A','SD2112',
             'TD1','TMW1_112','TMW1_656','lpuph','mlc3','LR1',
             'LR10','LR11','LR12','LR13','LR14',
             'LR17','LR18','LR19','LR2','LR3',
             'LR4','LR6','LR7','LR8','LR9'
            ]

os.chdir('../../ComplementaryData/Step_02_DraftModels/other_Lreuteri_seq/')
for i in name_list:
    print('-----  draft models %s  -----'%i)
    faa_file = i + '.faa'
    comd_str = 'carve ' + faa_file+ ' -o ../other_Lreuteri_moldebuilding/CarveMe/'+ i+'.xml'
    os.system(comd_str)