#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Step_seq_RAST_info.py
:description : script
:param : 
:returns: 
:rtype: 
"""

import os

import pandas as pd
from Bio import SeqIO

os.chdir('../../ComplementaryData/Step_01_Sequences_analysis/')

rast_functions = 'Lreuteri_biogaia_v03/RAST/RAST_functions.tsv'
rast_subsystem = 'Lreuteri_biogaia_v03/RAST/RAST_subsystem.tsv'

proka_seq = 'Lreuteri_biogaia_v03/Lreuteri_biogaia_v03.gbk'

rast_functions_df = pd.read_csv(rast_functions, sep='\t', )
rast_subsystem_df = pd.read_csv(rast_subsystem, sep='\t')

COG_functions = '../Initial_data/cog-20.cog.csv'
COG_classes = '../Initial_data/cog-20.def.tab'

prokka_functions = 'Lreuteri_biogaia_v03/prokka_functions.tsv'

# %% prokka functions
id_set = set()
locus_tag = 'locus_tag'
gene_functions = open(prokka_functions, "w")
gene_functions.write("%s\t%s\t%s\t%s\t%s\t\n" %
                     ('Locus_tag',
                      'Protein_id',
                      'Start',
                      'Stop',
                      'Product'))

for seq in SeqIO.parse(proka_seq, "genbank"):
    for seq_feature in seq.features:
        geneSeq = seq_feature.extract(seq.seq)
        if seq_feature.type == "CDS":
            write_seq = True
            try:
                if seq_feature.qualifiers['locus_tag'][0] in id_set:
                    write_seq = False
                    print('%s\t\tduplicate in seq, locus_tag is %s' % (
                        seq_feature.qualifiers[locus_tag][0], seq_feature.qualifiers['locus_tag'][0]))
            except KeyError:
                print('%s\t\t have no %s' % (seq_feature.qualifiers['locus_tag'][0], locus_tag))
                continue

            gene_functions.write(seq_feature.qualifiers[locus_tag][0] + '\t')
            id_set.add(seq_feature.qualifiers[locus_tag][0])

            if 'inference' in seq_feature.qualifiers.keys():
                gene_functions.write(
                    seq_feature.qualifiers['inference'][-1].split(':')[-1] + '\t')

            gene_functions.write(str(seq_feature.location.start + 1) + '\t' + str(seq_feature.location.end) + '\t')

            gene_functions.write("%s\t\n" % (
                seq_feature.qualifiers['product'][0].split(',')[0],))
gene_functions.close()
prokka_df = pd.read_csv(prokka_functions, sep='\t', )

gene_functions_df = prokka_df

# %% merge COG functions

# COG_functions_df = pd.read_csv(COG_functions,sep = ',',header = None,names=[i for i in range(0, 13)])
COG_functions_new_file = 'Lreuteri_biogaia_v03/COG_functions.tsv'
COG_functions_new = open(COG_functions_new_file, "w")

protiens = set(gene_functions_df['Protein_id'].values)
for line in open(COG_functions, 'r'):
    protien_id = line.split(',')[2]
    if protien_id in protiens:
        COG_functions_new.write(line)

COG_functions_new.close()

# COG_functions_df = COG_functions_df[]
COG_functions_df = pd.read_csv(COG_functions_new_file, sep=',', header=None, names=[i for i in range(0, 13)])
COG_functions_df = COG_functions_df[[2, 6]]
COG_functions_df.columns = ['Protein_id', 'COG ID']
gene_functions_df = gene_functions_df.merge(COG_functions_df, how='left', on='Protein_id')

COG_classes_df = pd.read_csv(COG_classes, sep='\t', header=None, names=['COG ID',
                                                                        'COG functional category',
                                                                        'COG name', 'Gene associated with the COG',
                                                                        'Functional pathway associated with the COG',
                                                                        'PubMed ID',
                                                                        'PDB ID',
                                                                        ], encoding='cp1252')

gene_functions_df = gene_functions_df.merge(COG_classes_df, how='left', on='COG ID', )

# %% merge RAST functions

gene_functions_df = gene_functions_df.merge(rast_functions_df, how='left', on=['Start'])

# %%
import cobra

model = cobra.io.load_json_model('../../ModelFiles/iHL622.json')

# gene_in_model_df = pd.DataFrame([i.id for i in model.genes],columns=['Locus_tag'])
gene_list = [i.id for i in model.genes]
#
# gene_in_model_df = gene_in_model_df.merge(gene_functions_df,how = 'left',on =['Locus_tag'] )
# gene_in_model_df.to_csv('Lreuteri_biogaia_v03/gene_in_model_functions.tsv',sep = '\t')

gene_functions_df['In_model'] = gene_functions_df['Locus_tag'].isin(gene_list)

# %%

gene_functions_df.to_csv('Lreuteri_biogaia_v03/gene_functions.tsv', sep='\t')
