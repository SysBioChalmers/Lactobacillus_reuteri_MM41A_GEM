#!/usr/bin/env python
#-*- coding: utf-8 -*-
# Created by lhao at 2019-02-27


import os
from Bio.Blast import NCBIXML

os.chdir('/Users/lhao/Box Sync/Project_Lreuteri/02_DraftModels/Template/')

t_ids = ['iBT721','iNF517','iMP429','iYO844','iML1515']
t_models = ['models/'+i+'.xml'for i in t_ids]
t_seqs = ['seqs/'+i+'.faa'for i in t_ids]


os.system('makeblastdb -in Blast/final_6475_prot.faa -dbtype prot -out Blast/Lreu/Lreu -parse_seqids')
for index in range(len(t_ids)):
    t_seqsdb = 'Blast/' + t_ids[index] + '/' + t_ids[index]

    balstdb = 'makeblastdb -in ' + t_seqs[index] + ' -dbtype prot -out ' + t_seqsdb + ' -parse_seqids'
    #os.system(balstdb)

    outfile1 = 'Blast/'+ t_ids[index] +'_in_Lreu.xml'
    balsttoLeu = 'blastp -db Blast/Lreu/Lreu -query ' + t_seqs[index] + ' -out ' + outfile1 +  ' -outfmt 5 -evalue 1e-10 -max_target_seqs 1 '
    #os.system(balsttoLeu)


    outfile2 = 'Blast/Lreu_in_' + t_ids[index] + '.xml'
    balstfromLeu = 'blastp -db '+ t_seqsdb +' -query Blast/Final_6475_prot.faa -out ' + outfile2 + ' -outfmt 5 -evalue 1e-10 -max_target_seqs 1 '
    os.system(balstfromLeu)

#os.system('blastp -db Llac1363/Llac1363 -query Final_6475_prot.faa -out reu_in_lac.xml -outfmt 5 -evalue 1e-10 -max_target_seqs 1 ')
#os.system('blastp -db Lreu6475/Lreu6475 -query Llac1363_changedid_tag.faa -out lac_in_reu.xml -outfmt 5 -evalue 1e-10 -max_target_seqs 1 ')


#matched Genes:
# set cutoff
#   maxE              only look at genes with E-values <= this value (opt,
#                     default 10^-30)
#   minLen            only look at genes with alignment length >= this
#                     value (opt, default 200)
'question about minLen !!!'
#   minIde            only look at genes with identity >= this value
#                     (opt, default 40 (%))

'''
maxE = 10 ** -30
minLen = 200
minIde = 40


files_list =  ["reu_in_lac.xml","lac_in_reu.xml"]

set_ft = set()
set_tf = set()

for i in range(0,len(files_list)):
    result_handle = open(files_list[i])
    blast_records = NCBIXML.parse(result_handle)

    for blast_record in blast_records:
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                if hsp.expect <= 0.4 and hsp.identities >= 40 and hsp.align_length >= 200:
                    from_id =  blast_record.query.split(' ')[0]
                    to_id = alignment.hit_id
                    if i ==0:
                        set_ft.add(from_id+'\t'+to_id)
                    else:
                        set_tf.add(to_id+'\t'+from_id)

matchlist = set_ft & set_tf

targetgene_list = [j.split('\t')[1] for j in matchlist ]
'''



'''
# read TXT
matchset1 = set()
for line in open('reu_in_lac.txt'):

    if '#' not in line:
        recordi  = line.split('\t')
        if float(recordi[2]) > 40 and float(recordi[10])<(0.01):
            pair = recordi[0], recordi[1]
            matchset1.add(pair)


matchset2 = set()
for line in open('lac_in_reu.txt'):

    if '#' not in line:
        recordi  = line.split('\t')
        if float(recordi[2]) > 40 and float(recordi[10])<(0.01):
            pair = recordi[1],recordi[0]
            matchset2.add(pair)

matchset = matchset1 & matchset2

outfile  = open('matchgeneslist.txt','w')
matchset1 = list(matchset1)
matchset2 = list(matchset2)
matchset = list(matchset)


for i in range(0,max(len(matchset1),len(matchset2))):
    line = ['']*7
    if i <len(matchset1):
        line[0] = matchset1[i][0]
        line[1] = matchset1[i][1]

    if i < len(matchset2):
        line[2] = matchset2[i][0]
        line[3] = matchset2[i][1]

    if i < len(matchset):
        line[4] = matchset[i][0]
        line[5] = matchset[i][1]
    outfile.write('\t'.join(line)+'\n')

outfile.close()

'''