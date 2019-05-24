#!/usr/bin/env python
#-*- coding: utf-8 -*-
# Created by lhao at 2019-02-25
#
# Blast three version seqs.
# output the BBH number

import os
from Bio.Blast import NCBIXML
import itertools

os.chdir('../../ComplementaryData/Step_01_Sequences_analysis/')

names = ['Lreuteri_refseq_v01','Lreuteri_refseq_v02','Lreuteri_biogaia_v03']

seqs_ids = [ i + '/' + i + '.faa' for i in names]
seqs_dbs = [ 'blast/' + i + '/' + i  for i in names]


#make balst database
for index in range(len(seqs_ids)):
    balstdb = 'makeblastdb -in ' + seqs_ids[index] + ' -dbtype prot -out ' + seqs_dbs[index] + ' -parse_seqids'
    os.system(balstdb);

#blast pairwise

for index_pair in itertools.permutations(range(len(seqs_ids)),2):
    seq = seqs_ids[index_pair[0]]
    db = seqs_dbs[index_pair[1]]
    outfile = 'blast/'+ names[index_pair[0]] +'_in_' + names[index_pair[1]] + '.csv'

    balstcmd = 'blastp -db ' + db + ' -query ' + seq + ' -out ' + outfile +  ' -evalue 10e-5  -outfmt "6 qseqid sseqid evalue pident length bitscore ppos"'
    os.system(balstcmd);


print(os.system('ls'))

