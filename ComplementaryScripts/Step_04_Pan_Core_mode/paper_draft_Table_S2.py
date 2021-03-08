#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 3/3/21

"""Table_S2.py
:description : script
:param : 
:returns: 
:rtype: 
"""

import os
import pandas as pd

from Bio import SeqIO


print('----- loading data -----')
os.chdir('../../ComplementaryData/Step_02_DraftModels/other_Lreuteri_seq/')

Lreu_seq = '../Lreuteri_biogaia_v03_2.faa'

name_list = ['100_23','20_2','3c6','ATCC53608','CF48_3A',
             'DSM200016','I5007','IRT','JCM1112',
             'LTH2584','LTH5448','MM2_3','MM4_1A','SD2112',
             'TD1','TMW1_112','TMW1_656','lpuph','mlc3','LR1',
             'LR10','LR11','LR12','LR13','LR14',
             'LR17','LR18','LR19','LR2','LR3',
             'LR4','LR6','LR7','LR8','LR9'
            ]
df = pd.DataFrame(index = name_list,columns=[
    'id','date','taxonomy','references','dbxrefs','annotations'])


for i in name_list:
    gbk_file = i + '.gbff'


    # record = SeqIO.read(gbk_file, "genbank")
    for record in SeqIO.parse(gbk_file, "genbank"):
        record.id
        record.name
        record.letter_annotations
        record.dbxrefs
        record.annotations["source"]


        df.loc[i] = [record.id,record.annotations['date'],record.annotations['taxonomy'],
                     record.annotations['references'][0].title,
                     record.dbxrefs,record.annotations,]
        break


df.to_csv('table_s1_.tsv',sep = '\t')