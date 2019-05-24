#!/usr/bin/env python
#-*- coding: utf-8 -*-
# Created by lhao at 2019-05-16

import os
import itertools
import My_def
from matplotlib_venn import venn3, venn3_circles
from matplotlib import pyplot as plt

os.chdir('../../ComplementaryData/Step_01_Sequences_analysis/')

names = ['Lreuteri_refseq_v01','Lreuteri_refseq_v02','Lreuteri_biogaia_v03']
result_dic = {}
for index_pair in itertools.combinations(range(len(names)),2):
    result1 = 'blast/' + names[index_pair[0]] + '_in_' + names[index_pair[1]] + '.csv'
    result2 = 'blast/' + names[index_pair[1]] + '_in_' + names[index_pair[0]] + '.csv'
    resultdf = My_def.select_blast(result1, result2, best_match=True,evalue = 10**-10, pident = 0, length = 0, bitscore = 100, ppos = 45)
    result_dic[str(index_pair[0]+1) + ' and ' + str(index_pair[1]+1) + ' BBH'  ] = len(resultdf)

print(result_dic)






#TODO：fig？
'''
textlist = []
v = venn3(subsets=(1, 1, 1, 1, 1, 1, 0))
listid = ['100', '010', '110', '001', '101', '011', '111']
for i in range(7):
    v.get_label_by_id(listid[i]).set_text(textlist[i])

'''




