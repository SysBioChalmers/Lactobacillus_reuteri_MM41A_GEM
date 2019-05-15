#!/usr/bin/env python
#-*- coding: utf-8 -*-
# Created by lhao at 2019-02-26
# useless!!!
import os

os.chdir('/Users/lhao/Box Sync/Project_Lreuteri/02_DraftModels/Template/')

uniportfile_name  = 'uniprot-Lactococcus+lactis+subsp.+cremoris+(strain+MG1363).txt'

outlist_file = open('Llca1363_mapgenelist.txt','w')

for line in open(uniportfile_name):
    line = line.replace('\n','')

    if line.startswith('ID   '):
        #find gene_i id
        gene_i_list = ['']*10
        temp = line.split('ID   ')[1]
        gene_i_id = temp.split(' ')[0]
        gene_i_list[0] = gene_i_id
    elif line.startswith('AC   '):
        temp = line.split('AC   ')[1]
        gene_i_AC = temp.split(';')[0]
        gene_i_list[1] = gene_i_AC
    elif line.startswith('DR   RefSeq; '):
        temp = line.split('DR   RefSeq; ')[1]
        temp = temp.replace(' ','')

        gene_i_refid = temp
        if ';' in temp:
            gene_i_refid_simple = gene_i_refid.split(';')[0]
        else :
            gene_i_refid_simple = gene_i_refid

        gene_i_list[2] = gene_i_refid_simple
        gene_i_list[3] = gene_i_refid


    elif line.startswith('DR   PATRIC; '):
        temp = line.split('DR   PATRIC; ')[1]
        temp = temp.replace(' ', '')
        if ';' in temp:
            gene_i_patric = temp.split(';')[0]
        else:
            gene_i_patric = temp
        gene_i_patric_simple = 'peg' + gene_i_patric.split('peg.')[1]
        gene_i_list[4] = gene_i_patric_simple
        gene_i_list[5] = gene_i_patric
    elif line.startswith('DR   KEGG; '):
        temp = line.split('DR   KEGG; ')[1]
        temp = temp.replace(' ', '')
        if ';' in temp:
            gene_i_kegg = temp.split(';')[0]
        else:
            gene_i_kegg = temp

        gene_i_kegg_simple = gene_i_kegg.split(":")[1]
        gene_i_list[6] = gene_i_kegg_simple
        gene_i_list[7] = gene_i_kegg

    elif line.startswith('//'):

        outlist_file.write('\t'.join(gene_i_list)+'\t\n')

outlist_file.close()



