#!/usr/bin/env python
#-*- coding: utf-8 -*-
# Created by lhao at 2019-02-27


import os
import My_def

os.chdir('../../ComplementaryData/02_DraftModels/Template/')

t_ids = ['iBT721','iNF517','iMP429','iYO844','iML1515']
t_models = ['template_models/'+i+'.xml'for i in t_ids]


Lreu_seq = '../Lreuteri_biogaia_v03.faa'
Lreu_db = 'blast/Lreu/Lreu'
os.system('makeblastdb -in ../Lreuteri_biogaia_v03.faa -dbtype prot -out blast/Lreu/Lreu -parse_seqids')

for index in range(len(t_ids)):  #
    t_seqsdb = 'blast/' + t_ids[index] + '/' + t_ids[index]
    t_seqs = 'template_seqs/' + t_ids[index]+'.faa'

    balstdb = 'makeblastdb -in ' + t_seqs + ' -dbtype prot -out ' + t_seqsdb + ' -parse_seqids'
    os.system(balstdb)

    outfile1 = 'blast/'+ t_ids[index] +'_in_Lreu.csv'
    balsttoLeu = 'blastp -db Blast/Lreu/Lreu -query ' + t_seqs + ' -out ' + outfile1 +  ' -evalue 10e-5  -outfmt "6 qseqid sseqid evalue pident length bitscore ppos"'
    os.system(balsttoLeu)

    outfile2 = 'blast/Lreu_in_' + t_ids[index] + '.csv'
    balstfromLeu = 'blastp -db '+ t_seqsdb +' -query ../Lreuteri_biogaia_v03.faa -out ' + outfile2 + '  -evalue 10e-5  -outfmt "6 qseqid sseqid evalue pident length bitscore ppos"'
    os.system(balstfromLeu)

    result_df = My_def.select_blast(outfile1, outfile2, best_match=True,evalue = 10**-10, pident = 0, length = 0, bitscore = 100, ppos = 45)

    result_df.to_csv('blast/' + t_ids[index] + '_and_Lreu.csv' , index = False)


