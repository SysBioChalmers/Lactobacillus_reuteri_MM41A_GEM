#!/usr/bin/env python
#-*- coding: utf-8 -*-
# Created by lhao at 2019-02-27


import os
from Bio import SeqIO


os.chdir('/Users/lhao/Box Sync/Project_Lreuteri/02_DraftModels/Template/seqs/')

templatelist = ['iBT721',
                'iNF517',
                'iMP429',
                'iYO844',
                'iML1515']
templatelist = ['iML1515'
                ]

for i in templatelist:
    gbk_file = i + '.gbff'
    fna_file = i + '.fna'
    faa_file = i + '.faa'
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
        locus_tag = 'locus_tag'
        if i in ['iMP429', 'iYO844']:
            locus_tag = 'old_locus_tag'

        for seq_feature in seq.features :
            geneSeq = seq_feature.extract(seq.seq)
            if seq_feature.type == "CDS" :
                try:

                    geneNC.write(">%s %s, %s\n%s\n" % (
                        seq_feature.qualifiers[locus_tag][0],
                        seq_feature.qualifiers['protein_id'][0],
                        seq.description,
                        geneSeq))

                    assert len(seq_feature.qualifiers['translation']) == 1
                    geneAA.write(">%s %s, %s\n%s\n" % (
                        seq_feature.qualifiers[locus_tag][0],
                        seq_feature.qualifiers['protein_id'][0],
                        seq.description,
                        seq_feature.qualifiers['translation'][0]))
                    listfile.write("%s\t%s\t\t\n" % (
                        seq_feature.qualifiers[locus_tag][0],
                        seq_feature.qualifiers['protein_id'][0],))

                except:
                    print(seq_feature.qualifiers['locus_tag'][0])

    input.close()
    geneNC.close()
    geneAA.close()
    listfile.close()



