#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 2/12/20

"""Step_04_1_get_seq_files_for_pan.py
:description : script
:param : 
:returns: 
:rtype: 
"""

import os

import pandas as pd

os.chdir('../../ComplementaryData/Initial_data/other_Lreuteri_seq/')
# %% <get seqs from PATRIC>
# cmd = 'source "/Applications/PATRIC.app//user-env.zsh" \n ' \
#     'p3-all-genomes --eq "genome_name,Lactobacillus reuteri" | p3-get-genome-data >Lreuteri_summery.tbl'
# os.system(cmd)

Lreuteri_summery = pd.read_csv('Lreuteri_summery.tbl', sep='\t', header=0, dtype='str')
genome_ids = Lreuteri_summery['genome.genome_id'].to_list()

ftp_path = 'ftp://ftp.patricbrc.org/genomes_by_species/Lactobacillus_reuteri/'

file_types = ['.fna', '.PATRIC.faa', '.PATRIC.features.tab', '.PATRIC.subsystem.tab', '.PATRIC.pathway.tab',
              '.RefSeq.faa', '.RefSeq.features.tab', ]
# errors = open('download_errors.txt','w')
all_seq_set = set(os.listdir('seqs_from_PATRIC/'))

for index in range(0, len(genome_ids)):
    id = genome_ids[index]
    for file_type in file_types:
        try:
            if str(id) + file_type not in all_seq_set:
                file_url = ftp_path + str(id) + '/' + str(id) + file_type
                # 'ftp://ftp.patricbrc.org/genomes_by_species/Lactobacillus_reuteri//349123.13/349123.13.RedSeq.faa'

                # wget.download(file_url, 'seqs_from_PATRIC/')
        except:
            # raise
            print('ERROR!!!   ' + file_url)
            # errors.write('ERROR!!!   ' + file_url +'\n')
    print(index)
all_seq_set = set(os.listdir('seqs_from_PATRIC/'))
ref_seq_set = [i for i in all_seq_set if 'Ref' in i]

with_file = [*map(lambda x: x.split('.')[0] + '.' + x.split('.')[1], all_seq_set)]
with_reference = [*map(lambda x: x.split('.')[0] + '.' + x.split('.')[1], ref_seq_set)]

# Lreuteri_summery['genome.genome_id.1'].to_string()
# Lreuteri_summery.loc[Lreuteri_summery['genome.genome_id'].isin(with_file),'seq_files'] = 'yes'
Lreuteri_summery['seq_files'] = Lreuteri_summery['genome.genome_id'].isin(with_file)
Lreuteri_summery['Refseq_files'] = Lreuteri_summery['genome.genome_id'].isin(with_reference)

Lreuteri_summery['Strain_name'] = Lreuteri_summery['genome.genome_name'].str.replace('Lactobacillus reuteri ', '')
Lreuteri_summery['Strain_name'] = Lreuteri_summery['Strain_name'].str.replace('strain ', '')
Lreuteri_summery['Strain_name'] = Lreuteri_summery['Strain_name'].str.replace(' ', '_')
Lreuteri_summery = Lreuteri_summery.sort_values(by=['seq_files', 'Refseq_files', 'Strain_name'], ascending=False)
Lreuteri_summery.to_csv('Lreuteri_summery.csv', sep='\t', index=False)

# %%  <get seqs from NCBI>
# refseq_database = pd.read_csv('assembly_summary_refseq.txt',sep = '\t',header=1)
# refseq_Lreu = refseq_database[refseq_database.organism_name.str.contains('Lactobacillus reuteri')]
#
# refseq_Lreu.to_csv('Lreu_assembly_summary_refseq.txt',sep = '\t',index=False)

# refseq_Lreu.columns
# refseq_Lreu['# assembly_accession']
# refseq_Lreu.ftp_path

# %%
#
# for index in range(0,refseq_Lreu.shape[0]):
#     ftp_path = refseq_Lreu.iloc[index]['ftp_path']
#     file_url = ftp_path  + '/' + ftp_path.split('/')[-1] + '_genomic.gbff.gz'
#     # strain_name = refseq_Lreu.iloc[index]['organism_name']
#     # if strain_name == 'Lactobacillus reuteri':
#     #     strain_name = 'Lactobacillus reuteri ' + refseq_Lreu.iloc[index]['infraspecific_name'].split('=')[-1]
#     strain_name = refseq_Lreu.iloc[index]['infraspecific_name'].split('=')[-1]
#     try:
#         # wget.download(file_url,'all_Lreu_strains_gbff/' )
#         # os.rename('all_Lreu_strains_gbff/'+ftp_path.split('/')[-1] + '_genomic.gbff.gz', os.path.join('all_Lreu_strains_gbff/', strain_name+ '.gbff.gz'))
#
#         os.system('gzip -d ' +  'seqs_from_NCBI/'+strain_name.replace(' ','\ ') + '.gbff.gz')
#         # os.rename('all_Lreu_strains_gbff/'+ftp_path.split('/')[-1] + '_genomic.gbff', os.path.join('all_Lreu_strains_gbff/', strain_name+ '.gbff'))
#         # os.remove(my_file)
#         print(index)
#     except:
#         raise
#         print('ERROR!!!   ' + file_url)

# NOTE: I49 and DSM 20016 have two versions
