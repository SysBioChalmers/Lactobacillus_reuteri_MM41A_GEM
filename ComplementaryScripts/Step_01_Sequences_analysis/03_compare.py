#!/usr/bin/env python
#-*- coding: utf-8 -*-
# Created by lhao at 2019-05-16

import os
import itertools
import My_def
from matplotlib_venn import venn3, venn3_circles
from matplotlib import pyplot as plt
import pandas as pd

os.chdir('../../ComplementaryData/Step_01_Sequences_analysis/blast/')

#names = ['Lreuteri_refseq_v01','Lreuteri_refseq_v02','Lreuteri_biogaia_v03']
result12 = pd.read_csv('v01andv02.csv',usecols = [0,1])
result13 = pd.read_csv('v01andv03.csv',usecols = [0,1])
result23 = pd.read_csv('v02andv03.csv',usecols = [0,1])

result12.columns = ['v02', 'v01']
result13.columns = ['v03', 'v01']
result23.columns = ['v03', 'v02']



# %% option1
result = pd.merge(result13,result23,on=['v03'],how='inner')
result = pd.merge(result,result12,on=['v02','v01'],how='inner')
result = result.fillna(0)



a111 = len(result)

a110 = len(result12) - a111
a101 = len(result13) - a111
a011 = len(result23) - a111

a100 = 2019- a111 - a110 - a101
a010 = 1919 -a111 - a110 - a011
a001 = 2019 -a111 - a101 - a011





# %% option2 (filed)
# result = pd.merge(result13,result23,on=['v03'],how='outer')
# result = pd.merge(result,result12,on=['v02','v01'],how='outer')
# result = result.fillna(0)
#
#
# a111 = len(result[(result['v01']!=0) & (result['v01']!=0) & (result['v01']!=0) ])
#
# a110 = len(result[(result['v01']!=0) & (result['v02']!=0) & (result['v03']==0) ])
# a101 = len(result[(result['v01']!=0) & (result['v02']==0) & (result['v03']==0) ])
# a011 = len(result[(result['v01']==0) & (result['v02']!=0) & (result['v03']==0) ])
#
# a100 = 2019- a111 - a110 - a101
# a010 = 1919 -a111 - a110 - a011
# a001 = 2019 -a111 - a101 - a011

# %% fig

textlist = (a100, a010, a110, a001, a101, a011, a111)
listid = ['100', '010', '110', '001', '101', '011', '111']

plt.figure(figsize=(4,4))

v = venn3(subsets=(1, 1, 1, 1, 1, 1, 1), set_labels = ('V1', 'V2', 'V3'))

for i in range(7):
    v.get_label_by_id(listid[i]).set_text(textlist[i])
    v.get_label_by_id(listid[i]).set_fontsize(14)

#for text in v.subset_labels:
#    text.set_fontsize(16)
plt.show()



