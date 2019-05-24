#!/usr/bin/env python
#-*- coding: utf-8 -*-
# Created by lhao at 2019-02-25
#
# annote the raw genmoe by Prokka

import os

os.chdir('../../ComplementaryData/Step_01_Sequences_analysis/')
#prokka without any protien sequences
#os.system('prokka --outdir prokka/v1 --prefix Lreuteri_6475_v1 Final_6475_Sequence.fas --force')

# trusted proteins from refseq
    #results: the same with above but protien name . better
os.system('prokka '
          ' --proteins Lreuteri_refseq_v02/Lreuteri_refseq_v02.faa'
          ' --outdir Lreuteri_biogaia_v03/prokka/ --prefix Lreuteri_biogaia_v03 Lreuteri_biogaia_v03/Lreuteri_biogaia_v03_Complete_sequence.fas'
          ' --force')

os.system('cp Lreuteri_biogaia_v03/prokka/Lreuteri_biogaia_v03.faa Lreuteri_biogaia_v03/Lreuteri_biogaia_v03.faa')
os.system('cp Lreuteri_biogaia_v03/prokka/Lreuteri_biogaia_v03.gbk Lreuteri_biogaia_v03/Lreuteri_biogaia_v03.gbk')

print(os.system('ls'))

