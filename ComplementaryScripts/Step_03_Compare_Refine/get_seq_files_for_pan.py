#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 2/12/20

"""get_seq_files_for_pan.py
:description : script
:param : 
:returns: 
:rtype: 
"""

import os
import pandas as pd
import wget
import gzip

os.chdir('../../ComplementaryData/Initial_data/other_Lreuteri_seq/')

refseq_database = pd.read_csv('assembly_summary_refseq.txt',sep = '\t',header=1)
refseq_Lreu = refseq_database[refseq_database.organism_name.str.contains('Lactobacillus reuteri')]

refseq_Lreu.to_csv('Lreu_assembly_summary_refseq.txt',sep = '\t',index=False)

# refseq_Lreu.columns
# refseq_Lreu['# assembly_accession']
# refseq_Lreu.ftp_path

# %%

for index in range(0,refseq_Lreu.shape[0]):
    ftp_path = refseq_Lreu.iloc[index]['ftp_path']
    file_url = ftp_path  + '/' + ftp_path.split('/')[-1] + '_genomic.gbff.gz'
    # strain_name = refseq_Lreu.iloc[index]['organism_name']
    # if strain_name == 'Lactobacillus reuteri':
    #     strain_name = 'Lactobacillus reuteri ' + refseq_Lreu.iloc[index]['infraspecific_name'].split('=')[-1]
    strain_name = refseq_Lreu.iloc[index]['infraspecific_name'].split('=')[-1]
    try:
        # wget.download(file_url,'all_Lreu_strains_gbff/' )
        # os.rename('all_Lreu_strains_gbff/'+ftp_path.split('/')[-1] + '_genomic.gbff.gz', os.path.join('all_Lreu_strains_gbff/', strain_name+ '.gbff.gz'))

        os.system('gzip -d ' +  'all_Lreu_strains_gbff/'+strain_name.replace(' ','\ ') + '.gbff.gz')
        # os.rename('all_Lreu_strains_gbff/'+ftp_path.split('/')[-1] + '_genomic.gbff', os.path.join('all_Lreu_strains_gbff/', strain_name+ '.gbff'))
        # os.remove(my_file)
        print(index)
    except:
        raise
        print('ERROR!!!   ' + file_url)

# NOTE: I49 and DSM 20016 have two versions

