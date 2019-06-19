#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 2019-06-18

"""seq_analysis.py
:description : script to process the seqs of Lreuteri
:param : input three versions of genmoes
:returns: out put processed seqs
:rtype: 
"""

import os
from Bio import SeqIO
import My_def
import pandas as pd

# %% <step: standlize the v01 and v02 version Lreu >

# move file from iterinal dir to Step_01_Sequences_analysis
initial_dir = 'ComplementaryData/Initial_data/Lreuteri_genomes/'
os.chdir(initial_dir)
'GCF_000159475.2_ASM15947v2_genomic.gbf'
os.system('cat Lreuteri_refseq_v01/NZ_ACGX00000000.scaffold.gbk/*.gbk >'+'../../Step_01_Sequences_analysis/Lreuteri_refseq_v01/'+'Lreuteri_refseq_v01.gbk')
os.system('cp Lreuteri_refseq_v02/GCF_000159475.2_ASM15947v2_genomic.gbff ../../Step_01_Sequences_analysis/Lreuteri_refseq_v02/Lreuteri_refseq_v02.gbff')
os.system('cp Lreuteri_biogaia_v03.fas ../../Step_01_Sequences_analysis/Lreuteri_biogaia_v03/Lreuteri_biogaia_v03.fas')

os.chdir('../../Step_01_Sequences_analysis/')

# standlize and remove duplications
My_def.seq_ana.gbk2faa('Lreuteri_refseq_v01/Lreuteri_refseq_v01.gbk','Lreuteri_refseq_v01/Lreuteri_refseq_v01.faa','protein_id')
My_def.seq_ana.gbk2faa('Lreuteri_refseq_v02/Lreuteri_refseq_v02.gbff','Lreuteri_refseq_v02/Lreuteri_refseq_v02.faa','protein_id')

os.chdir('../../')


# %% <step: annotion the v03 seq by prokka>
os.chdir('ComplementaryData/Step_01_Sequences_analysis/')
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
os.chdir('../../')

# %% <step: blast_pairwise>
os.chdir('ComplementaryData/Step_01_Sequences_analysis/')

seq1 = 'Lreuteri_refseq_v01/Lreuteri_refseq_v01.faa'
seq2 = 'Lreuteri_refseq_v02/Lreuteri_refseq_v02.faa'
seq3 = 'Lreuteri_biogaia_v03/Lreuteri_biogaia_v03.faa'


qseq_file = seq1
sseq_file = seq2

blast_result_s_in_q, blast_result_q_in_s = My_def.seq_ana.blastp_pairwise(qseq_file, sseq_file, out_dir='blast/',diamond=True)
blast_result_df12 = My_def.seq_ana.select_blast(blast_result_s_in_q, blast_result_q_in_s, best_match=True, evalue=10 ** -10, pident=0, length=0,
                                      bitscore=100, ppos=45, qcovs=0)

blast_result_df12.to_csv('blast/' + 'v01andv02.csv', index=False)
os.system('rm ' + blast_result_s_in_q )
os.system('rm ' + blast_result_q_in_s)


qseq_file = seq1
sseq_file = seq3

blast_result_s_in_q, blast_result_q_in_s = My_def.seq_ana.blastp_pairwise(qseq_file, sseq_file, out_dir='',diamond=True)
blast_result_df13 = My_def.seq_ana.select_blast(blast_result_s_in_q, blast_result_q_in_s, best_match=True, evalue=10 ** -10,
                                      pident=0, length=0,
                                      bitscore=100, ppos=45, qcovs=0)

blast_result_df13.to_csv('blast/' + 'v01andv03.csv', index=False)
os.system('rm ' + blast_result_s_in_q)
os.system('rm ' + blast_result_q_in_s)


qseq_file = seq2
sseq_file = seq3

blast_result_s_in_q, blast_result_q_in_s = My_def.seq_ana.blastp_pairwise(qseq_file, sseq_file, out_dir='',diamond=True)
blast_result_df23 = My_def.seq_ana.select_blast(blast_result_s_in_q, blast_result_q_in_s, best_match=True, evalue=10 ** -10,
                                      pident=0, length=0,
                                      bitscore=100, ppos=45, qcovs=0)

blast_result_df23.to_csv('blast/' + 'v02andv03.csv', index=False)
os.system('rm ' + blast_result_s_in_q)
os.system('rm ' + blast_result_q_in_s)

# %% <step: comapration and print a venn fig>

from matplotlib_venn import venn3, venn3_circles
from matplotlib import pyplot as plt


seq1 = 'Lreuteri_refseq_v01/Lreuteri_refseq_v01.faa'
seq2 = 'Lreuteri_refseq_v02/Lreuteri_refseq_v02.faa'
seq3 = 'Lreuteri_biogaia_v03/Lreuteri_biogaia_v03.faa'
un1 = len(SeqIO.index(seq1, "fasta"))
un2 = len(SeqIO.index(seq2, "fasta"))
un3 = len(SeqIO.index(seq3, "fasta"))

result12 = pd.read_csv('blast/v01andv02.csv',usecols = [0,1])
result13 = pd.read_csv('blast/v01andv03.csv',usecols = [0,1])
result23 = pd.read_csv('blast/v02andv03.csv',usecols = [0,1])

result12.columns = ['v02', 'v01']
result13.columns = ['v03', 'v01']
result23.columns = ['v03', 'v02']



# %% <option1 >

result = pd.merge(result13,result23,on=['v03'],how='inner')
result = pd.merge(result,result12,on=['v02','v01'],how='inner')
result = result.fillna(0)

a111 = len(result)
a110 = len(result12) - a111
a101 = len(result13) - a111
a011 = len(result23) - a111
a100 = un1- a111 - a110 - a101
a010 = un2 -a111 - a110 - a011
a001 = un3 -a111 - a101 - a011

# %% option2 (filed)
# result = pd.merge(result13,result23,on=['v03'],how='outer')
# result = pd.merge(result,result12,on=['v02','v01'],how='outer')
# result = result.fillna(0)
#
#
# a111 = len(result[(result['Lreuteri_refseq_v01']!=0) & (result['Lreuteri_refseq_v01']!=0) & (result['Lreuteri_refseq_v01']!=0) ])
#
# a110 = len(result[(result['Lreuteri_refseq_v01']!=0) & (result['Lreuteri_refseq_v02']!=0) & (result['v03']==0) ])
# a101 = len(result[(result['Lreuteri_refseq_v01']!=0) & (result['Lreuteri_refseq_v02']==0) & (result['v03']==0) ])
# a011 = len(result[(result['Lreuteri_refseq_v01']==0) & (result['Lreuteri_refseq_v02']!=0) & (result['v03']==0) ])
#
# a100 = 2019- a111 - a110 - a101
# a010 = 1919 -a111 - a110 - a011
# a001 = 2019 -a111 - a101 - a011

# %% fig

textlist = (a100, a010, a110, a001, a101, a011, a111)
listid = ['100', '010', '110', '001', '101', '011', '111']

plt.figure(figsize=(4,4))

v = venn3(subsets=(1, 1, 1, 1, 1, 1, 1), set_labels = ('V1', 'V2', 'V3'))

for i in range(7):
    v.get_label_by_id(listid[i]).set_text(textlist[i])
    v.get_label_by_id(listid[i]).set_fontsize(14)

#for text in v.subset_labels:
#    text.set_fontsize(16)
plt.savefig('Venn_of_Lreu_seqs.png')
plt.show()

#TODO: result is different when using dinomd to blast



