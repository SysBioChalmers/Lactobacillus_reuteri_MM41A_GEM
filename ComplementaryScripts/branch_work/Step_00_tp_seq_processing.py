#!/usr/bin/env python
#-*- coding: utf-8 -*-
# Created by lhao at 2019-02-27

#change the seq reads id

import os
import My_def
os.chdir('../../ComplementaryData/Step_02_DraftModels/Template/template_seqs/')

# %%
os.system('cp -r ../../../Initial_data/template_seqs/ ./')

templatelist = ['iBT721',
                'iNF517',
                'iMP429',
                'iYO844',
                'iML1515']

for i in templatelist:
    gbk_file = i + '.gbff'
    faa_file = i + '.faa'
    My_def.seq_ana.gbk2faa(gbk_file,faa_file,locus_tag = 'locus_tag')
    

