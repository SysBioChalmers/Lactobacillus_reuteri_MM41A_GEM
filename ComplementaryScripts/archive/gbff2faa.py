#!/usr/bin/env python
#-*- coding: utf-8 -*-
# Created by lhao at 2019-02-27


import os
from Bio import SeqIO

os.chdir('/Users/lhao/Box Sync/Project_Lreuteri/Step_02_DraftModels/Template/')


gbk_file = "GCF_000009425.1_ASM942v1_genomic.gbff"
fna_file = "Llac1363_changedid_tag.fna"
faa_file = "Llac1363_changedid_tag.faa"
test_file = "list.txt"
input = open(gbk_file, "r")
geneNC = open(fna_file, "w")
geneAA = open(faa_file, "w")
listfile = open(test_file,'w')

for seq in SeqIO.parse(input, "genbank") :
    print ("Dealing with GenBank file of %s, \nOutput: \ngene nc: %s \ngene aa: %s" % (
        seq.id,
        fna_file,
        faa_file))

    for seq_feature in seq.features :
        geneSeq = seq_feature.extract(seq.seq)
        if seq_feature.type == "CDS" :
            try:
                geneNC.write(">%s %s, %s\n%s\n" % (
                    seq_feature.qualifiers['locus_tag'][0],
                    seq_feature.qualifiers['protein_id'][0],
                    seq.description,
                    geneSeq))

                assert len(seq_feature.qualifiers['translation']) == 1
                geneAA.write(">%s %s, %s\n%s\n" % (
                    seq_feature.qualifiers['locus_tag'][0],
                    seq_feature.qualifiers['protein_id'][0],
                    seq.description,
                    seq_feature.qualifiers['translation'][0]))
                listfile.write("%s\t%s\t\t\n" % (
                    seq_feature.qualifiers['locus_tag'][0],
                    seq_feature.qualifiers['protein_id'][0],

                ))

            except:
                print(seq_feature.qualifiers['locus_tag'][0])

input.close()
geneNC.close()
geneAA.close()
listfile.close()



