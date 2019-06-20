#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 2019-06-19

"""seq_ana.py
:description : script
:param : 
:returns: 
:rtype: 
"""
import os
from Bio import SeqIO
import pandas as pd
import re

def gbk2faa(input,output,locus_tag = 'locus_tag',remove_duplicate = True):
    '''
    :param input: opend gbk seq file
    :param output: faa file name and path
    :param locus_tag: chose seq id defult is 'locus_tag' (''protein_id'')
    :param remove_duplicate: remove duplication or not
    :return: none
    '''

    geneAA = open(output, "w")
    geneNC = open(output.replace('.faa','.fna'), "w")
    id_set = set()
    for seq in SeqIO.parse(input, "genbank") :
        print ("Dealing with GenBank file of %s, \nOutput: \ngene nc: %s \ngene aa: %s" %
               (seq.id,
                output.replace('.faa','.fna'),
                output))

        for seq_feature in seq.features :
            geneSeq = seq_feature.extract(seq.seq)
            if seq_feature.type == "CDS" :

                write_seq = True
                if remove_duplicate:
                    try:
                        if seq_feature.qualifiers[locus_tag][0] in id_set:
                            write_seq = False
                            print('%s\t\tduplicate in seq, locus_tag is %s'%(seq_feature.qualifiers[locus_tag][0],seq_feature.qualifiers['locus_tag'][0]))
                    except KeyError:
                        print('%s\t\t have no %s'%(seq_feature.qualifiers['locus_tag'][0],locus_tag))
                        continue
                if write_seq:
                    try:
                        geneNC.write(">%s %s, %s\n%s\n" % (
                            seq_feature.qualifiers[locus_tag][0],
                            seq_feature.qualifiers['product'][0],
                            seq.description,
                            geneSeq))

                        assert len(seq_feature.qualifiers['translation']) == 1
                        geneAA.write(">%s %s, %s\n%s\n" % (
                            seq_feature.qualifiers[locus_tag][0],
                            seq_feature.qualifiers['product'][0],
                            seq.description,
                            seq_feature.qualifiers['translation'][0]))
                        id_set.add(seq_feature.qualifiers[locus_tag][0])

                    except:
                        print('%s\t\tError processing'%seq_feature.qualifiers['locus_tag'][0])
    geneNC.close()
    geneAA.close()


def blastp_pairwise(qseq_file,sseq_file,out_dir = '',diamond = True):

    '''

    :param qseq_file:
    :param sseq_file:
    :param out_dir:
    :return:  two files
    '''

    os.system('mkdir blast_tmps')


    if diamond:
        mk_db = 'diamond makedb --in '+ qseq_file +' -d blast_tmps/qseq_db \n' \
                'diamond makedb --in '+ sseq_file +' -d blast_tmps/sseq_db '
        os.system(mk_db)
        print('diamond blasting...')
        options = ' --top 10 --more-sensitive '
        #options = ''
        diamond_blastp_cmd ='diamond blastp -d blast_tmps/qseq_db -q ' + sseq_file + options +' -o '+ out_dir +\
                            'blast_result_s_in_q.csv --evalue 1 --outfmt 6 qseqid sseqid evalue pident length bitscore ppos qcovhsp \n' \
                            'diamond blastp -d blast_tmps/sseq_db -q ' + qseq_file + options + ' -o '+ out_dir +\
                            'blast_result_q_in_s.csv --evalue 1 --outfmt 6 qseqid sseqid evalue pident length bitscore ppos qcovhsp'
        os.system(diamond_blastp_cmd)

    else:
        mk_db = 'makeblastdb -in '+ qseq_file +' -dbtype prot -out blast_tmps/qseq_db -parse_seqids\n' \
                'makeblastdb -in '+ sseq_file +' -dbtype prot -out blast_tmps/sseq_db -parse_seqids'
        os.system(mk_db)
        print('blasting...')
        blastp_cmd ='blastp -db blast_tmps/qseq_db -query ' + sseq_file + ' -out '+ out_dir +'blast_result_s_in_q.csv -evalue 1 -outfmt "6 qseqid sseqid evalue pident length bitscore ppos qcovs" \n' \
                    'blastp -db blast_tmps/sseq_db -query ' + qseq_file + ' -out '+ out_dir +'blast_result_q_in_s.csv -evalue 1 -outfmt "6 qseqid sseqid evalue pident length bitscore ppos qcovs"'
        os.system(blastp_cmd)

    os.system('rm -r blast_tmps')
    print('out put files are ' + out_dir + 'blast_result_s_in_q.csv and blast_result_q_in_s.csv')
    return out_dir + 'blast_result_s_in_q.csv', out_dir + 'blast_result_q_in_s.csv'


def select_blast(result1,result2,best_match=True,evalue = 10**-10, pident = 40, length = 200, bitscore = 0, ppos = 0 ,qcovs = 0):
    '''
    #    find BBH(Bidirectional Best Hits) from blast results
    #     balstcmd = 'blastp -db ' + db + ' -query ' + seq + ' -out ' + outfile +  ' -evalue 10e-5  -outfmt "6 qseqid sseqid evalue pident length bitscore ppos"'

    selsect the result hits from blast result
    :param result1: blast result file -outfmt "6 qseqid sseqid evalue pident length bitscore ppos"
    :param result2: last result file -outfmt "6 qseqid sseqid evalue pident length bitscore ppos"
    :param best_match: BBH(Bidirectional Best Hits) from blast results, best mached and  unique
    :param evalue:  evalue = 10**-10
    :param pident:  ident
    :param length:  matched length
    :param bitscore:
    :param ppos:
    :return:    a dataframe meet the params
    examlple:
            result1 = 'Lreuteri_refseq_v01_in_Lreuteri_refseq_v02.csv'
            result2 = 'Lreuteri_refseq_v02_in_Lreuteri_refseq_v01.csv'
            result_df = select_blast(result1, result2, best_match=True,evalue = 10**-10, pident = 40, length = 0, bitscore = 0, ppos = 0)
    '''
    names = ['qseqid', 'sseqid', 'evalue', 'pident', 'length', 'bitscore', 'ppos','qcovs']

    df1 = pd.read_csv(result1, sep = '\t', names = names)
    df2 = pd.read_csv(result2, sep = '\t', names = names)

    for index  in [0,1]:
        dfi = [df1,df2][index]
        dfi = dfi[(dfi["evalue"] <= evalue) & (dfi["pident"] >= pident) & (dfi["length"] >= length) & (dfi["bitscore"] >= bitscore) & (dfi["ppos"] >= ppos) & (dfi["qcovs"] >= qcovs)]
        dfi = dfi.copy()
        dfi['sseqid'] = dfi.sseqid.apply(lambda x: re.sub('(\|$)|(^.*?\|)','',x))
        dfi = dfi.sort_values(['qseqid', "evalue" , 'pident','bitscore'], ascending=[True,True,False,False])

        if best_match:
            dfi = dfi.drop_duplicates(subset='qseqid', keep='first')
        if index ==0:
            df1 = dfi
        else:
            df2 = dfi

    df2 = df2.rename(columns={'qseqid':'sseqid','sseqid':'qseqid'})

    result_df = pd.merge(df1, df2, on=['qseqid','sseqid'], how='inner')
    return result_df


if __name__ == '__main__':
    pass