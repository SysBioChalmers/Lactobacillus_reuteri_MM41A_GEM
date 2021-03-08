#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 2/1/21

"""MGS.py
:description : script to find MGS gene
:param : 
:returns: 
:rtype: 
"""


# GenBank: CAA0302616.1
import os
qseq_file = 'MGS.fasta'
# qseq_file = 'blast_temp/'+'100_23' + '.faa'
# sseq_file = ''
import My_def

os.chdir('../../ComplementaryData/Step_04_Pan_Core_model/')

name_list = ['100_23','20_2','3c6','ATCC53608','CF48_3A',
             'DSM200016','I5007','IRT','JCM1112',
             'LTH2584','LTH5448','MM2_3','MM4_1A','SD2112',
             'TD1','TMW1_112','TMW1_656','lpuph','mlc3','LR1',
             'LR10','LR11','LR12','LR13','LR14',
             'LR17','LR18','LR19','LR2','LR3',
             'LR4','LR6','LR7','LR8','LR9'
            ]

seq_dir = '../../ComplementaryData/Initial_data/other_Lreuteri_seq/archive/'


for i in name_list:
    os.system('mkdir blast_tmps')
    print('diamond blasting...')
    sseq_file = seq_dir+i + '.gbff'
    faa_file = 'blast_temp/'+i.replace('-','_') + '.faa'
    # My_def.seq_ana.gbk2faa(sseq_file, faa_file, locus_tag='locus_tag')

    mk_db = 'diamond makedb --in ' + faa_file + ' -d blast_tmps/sseq_db '
    os.system(mk_db)

    options = ' --top 1 --more-sensitive '
    #options = ''
    diamond_blastp_cmd ='diamond blastp -d blast_tmps/sseq_db -q ' + qseq_file + options +' -o '+ 'blast_temp/' +i+\
                        'blast_result_s_in_q.csv --evalue 1 --outfmt 6 qseqid sseqid evalue pident length bitscore ppos qcovhsp \n'
    os.system(diamond_blastp_cmd)