#!/usr/bin/env python
#-*- coding: utf-8 -*-
# Created by lhao at 2019-05-15

import os
from Bio import SeqIO

os.chdir('../../ComplementaryData/Step_01_Sequences_analysis/')

gbk_file = 'Lreuteri_refseq_v01/Lreuteri_refseq_v01.gbk'
fna_file = "Lreuteri_refseq_v01/Lreuteri_refseq_v01.fna"
faa_file = "Lreuteri_refseq_v01/Lreuteri_refseq_v01.faa"
input = open(gbk_file, "r")
geneNC = open(fna_file, "w")
geneAA = open(faa_file, "w")
id_set = set()

for seq in SeqIO.parse(input, "genbank") :
    print ("Dealing with GenBank file of %s, \nOutput: \ngene nc: %s \ngene aa: %s" % (
        seq.id,
        fna_file,
        faa_file))

    for seq_feature in seq.features :
        geneSeq = seq_feature.extract(seq.seq)
        if seq_feature.type == "CDS" :

            try:
                if seq_feature.qualifiers['protein_id'][0] not in id_set:
                    geneNC.write(">%s %s, %s\n%s\n" % (
                        seq_feature.qualifiers['protein_id'][0],
                        seq_feature.qualifiers['product'][0],
                        seq.description,
                        geneSeq))

                    assert len(seq_feature.qualifiers['translation']) == 1
                    geneAA.write(">%s %s, %s\n%s\n" % (
                        seq_feature.qualifiers['protein_id'][0],
                        seq_feature.qualifiers['product'][0],
                        seq.description,
                        seq_feature.qualifiers['translation'][0]))

                    id_set.add(seq_feature.qualifiers['protein_id'][0])

            except:
                print(seq_feature.qualifiers['locus_tag'][0])

input.close()
geneNC.close()
geneAA.close()
