#!/usr/bin/env python
#-*- coding: utf-8 -*-
# Created by lhao at 2019-02-25

import os

os.chdir('../../ComplementaryData/Step_01_Sequences_analysis/')

#prokka wiout any protien sequences
#os.system('prokka --outdir prokka/v1 --prefix Lreuteri_6475_v1 Final_6475_Sequence.fas --force')

# trusted proteins from refseq
    #result: same but proteins name
#os.system('prokka '
#          ' --proteins Lreuteri_6475_GCF_000159475.2_ASM15947v2/GCF_000159475.2_ASM15947v2_protein.faa'
#          ' --outdir prokka/v2 --prefix Lreuteri_6475_v2 Final_6475_Sequence.fas'
#          ' --force')

#os.system('cp prokka/v1/Lreuteri_6475_v1.faa Final_6475_prot.faa')

#os.system('carve Final_6475_prot.faa -o CarveMe/Lreuteri_6475_CarveMe_01.xml')


print(os.system('ls'))

