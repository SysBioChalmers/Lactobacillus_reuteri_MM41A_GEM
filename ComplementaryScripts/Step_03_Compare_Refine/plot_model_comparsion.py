#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 12/17/20

"""plot_model_comparsion.py
:description : script
:param : 
:returns: 
:rtype: 
"""
import os

import cobra
import matplotlib.pyplot as plt

os.chdir('../../ComplementaryData/Step_03_Compare_Refine/')
print('----- loading model -----')
iHL622 = cobra.io.load_json_model('../../ModelFiles/iHL622.json')
iNF517 = cobra.io.read_sbml_model('../Initial_data/template_models/iNF517.xml')
LbReueri = cobra.io.read_sbml_model('../Initial_data/template_models/Lreuteri_530.xml')
iBT721 = cobra.io.load_json_model('../Step_02_DraftModels/Template/template_models/iBT721_standlized.json')
# iML1515 = cobra.io.read_sbml_model('../Initial_data/template_models/iML1515.xml')

set1 = set([i.id for i in iHL622.reactions])
set2 = set([i.id for i in iNF517.reactions])
set3 = set([i.id for i in LbReueri.reactions])
set4 = set([i.id for i in iBT721.reactions])

a111 = len(set1 & set2 & set3)
a110 = len(set1 & set2 - set3)
a101 = len(set1 & set3 - set2)
a011 = len(set3 & set2 - set1)
a100 = len(set1 - set2 - set3)
a010 = len(set2 - set1 - set3)
a001 = len(set3 - set2 - set1)

textlist = (a100, a010, a110, a001, a101, a011, a111)
listid = ['100', '010', '110', '001', '101', '011', '111']
from matplotlib_venn import venn3

plt.figure(figsize=(4, 4))

v = venn3(subsets=(1, 1, 1, 1, 1, 1, 1), set_labels=('$i$HL622', '$i$NF517', 'LbReueri'))

for i in range(7):
    v.get_label_by_id(listid[i]).set_text(textlist[i])
    v.get_label_by_id(listid[i]).set_fontsize(14)

# for text in v.subset_labels:
#    text.set_fontsize(16)
# plt.savefig('Venn_of_Lreu_seqs.png')
plt.show()

f1 = open('iHL622_reactions.txt ', 'w')
f1.write('\n'.join(set1))
f1.close()

f1 = open('iNF517_reactions.txt ', 'w')
f1.write('\n'.join(set2))
f1.close()

f1 = open('LbReueri_reactions.txt ', 'w')
f1.write('\n'.join(set3))
f1.close()

f1 = open('iBT721_reactions.txt ', 'w')
f1.write('\n'.join(set4))
f1.close()

# web tool: http://bioinformatics.psb.ugent.be/webtools/Venn/
