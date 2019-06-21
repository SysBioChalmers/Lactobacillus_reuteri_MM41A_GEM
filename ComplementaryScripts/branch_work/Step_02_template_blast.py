#!/usr/bin/env python
#-*- coding: utf-8 -*-
# Created by lhao at 2019-02-27
#balset pairwise and get bbh file

import os
import pandas as pd
import re
import My_def



os.chdir('../../ComplementaryData/Step_02_DraftModels/Template/')

t_ids = ['iBT721','iNF517','iMP429','iYO844','iML1515']

for i in t_ids:  #

    qseq_file = '../Lreuteri_biogaia_v03.faa'
    sseq_file = 'template_seqs/' + i +'.faa'

    blast_result_s_in_q, blast_result_q_in_s = My_def.seq_ana.blastp_pairwise(qseq_file, sseq_file, out_dir='blast/')

    blast_result_df = My_def.select_blast(blast_result_s_in_q, blast_result_q_in_s, best_match=True, evalue=10 ** -10, pident=0, length=0,
                                          bitscore=100, ppos=45, qcovs=0)

    blast_result_df.to_csv('blast/' + i + '_and_Lreu.csv', index=False)
    os.system('rm ' + blast_result_s_in_q )
    os.system('rm ' + blast_result_q_in_s)


