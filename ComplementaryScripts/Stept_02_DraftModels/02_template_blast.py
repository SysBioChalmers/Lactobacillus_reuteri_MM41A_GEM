#!/usr/bin/env python
#-*- coding: utf-8 -*-
# Created by lhao at 2019-02-27


import os
import pandas as pd
import re
import My_def

def blastp_pairwise(qseq_file,sseq_file,out_dir = ''):

    '''

    :param qseq_file:
    :param sseq_file:
    :param out_dir:
    :return:  two files
    '''

    os.system('mkdir blast_tmps')

    mk_db = 'makeblastdb -in '+ qseq_file +' -dbtype prot -out blast_tmps/qseq_db -parse_seqids\n' \
            'makeblastdb -in '+ sseq_file +' -dbtype prot -out blast_tmps/sseq_db -parse_seqids'
    os.system(mk_db)
    print('blasting...')
    blastp_cmd ='blastp -db blast_tmps/qseq_db -query ' + sseq_file + ' -out '+ out_dir +' blast_result_s_in_q.csv -evalue 1 -outfmt "6 qseqid sseqid evalue pident length bitscore ppos qcovs" \n' \
                'blastp -db blast_tmps/sseq_db -query ' + qseq_file + ' -out '+ out_dir +' blast_result_q_in_s.csv -evalue 1 -outfmt "6 qseqid sseqid evalue pident length bitscore ppos qcovs"'
    os.system(blastp_cmd)

    #os.system('rm -r blast_tmps')
    return out_dir + 'blast_result_s_in_q.csv', out_dir + 'blast_result_q_in_s.csv'
    print('out put files are ' + out_dir + '/blast_result_s_in_q.csv and blast_result_q_in_s.csv')


def select_blast(result1 = 'blast_result_s_in_q.csv',result2 = 'blast_result_q_in_s.csv',best_match=True,evalue = 10**-10, pident = 40, length = 200, bitscore = 0, ppos = 0 ,qcovs = 0):
    '''
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
        dfi = dfi.sort_values(['qseqid', 'pident','bitscore'], ascending=[True,False,False])

        if best_match:
            dfi = dfi.drop_duplicates(subset='qseqid', keep='first')

        if index ==0:
            df1 = dfi
        else:
            df2 = dfi

    df2 = df2.rename(columns={'qseqid':'sseqid','sseqid':'qseqid'})

    result_df = pd.merge(df1, df2, on=['qseqid','sseqid'], how='inner')
    return result_df


if __name__ =='__main__':

    os.chdir('../../ComplementaryData/Stept_02_DraftModels/Template/')

    t_ids = ['iBT721','iNF517','iMP429','iYO844','iML1515']

    for i in t_ids:  #

        qseq_file = '../Lreuteri_biogaia_v03.faa'
        sseq_file = 'template_seqs/' + i +'.faa'

        blast_result_s_in_q, blast_result_q_in_s = blastp_pairwise(qseq_file, sseq_file, out_dir='')

        blast_result_df = My_def.select_blast(blast_result_s_in_q, blast_result_q_in_s, best_match=True, evalue=10 ** -10, pident=0, length=0,
                                              bitscore=100, ppos=45, qcovs=0)

        blast_result_df.to_csv('blast/' + i + '_and_Lreu.csv', index=False)
        os.system('rm ' + blast_result_s_in_q )
        os.system('rm ' + blast_result_q_in_s)


