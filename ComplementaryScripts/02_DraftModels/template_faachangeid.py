#!/usr/bin/env python
#-*- coding: utf-8 -*-
# Created by lhao at 2019-02-27
# ！！！useless

import os
from Bio import SeqIO

os.chdir('/Users/lhao/Box Sync/Project_Lreuteri/02_DraftModels/Template/')
out_file = open('Llac1363_change_id.faa','w')

old_loc = []
ref_id = []
for line in open('list.txt'):
    temp = line.split('\t')
    old_loc.append(temp[0])
    ref_id.append(temp[1])

for seq in SeqIO.parse('GCF_000009425.1_ASM942v1_protein.faa', "fasta") :
    if seq.id in ref_id:
        new_id = old_loc[ref_id.index(seq.id)]
    out_file.write('>'+new_id+'\n')
    out_file.write(str(seq.seq)+'\n')



out_file.close()










