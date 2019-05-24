#!/usr/bin/env python
#-*- coding: utf-8 -*-
# Created by lhao at 2019-05-15
#    find BBH(Bidirectional Best Hits) from blast results
#     balstcmd = 'blastp -db ' + db + ' -query ' + seq + ' -out ' + outfile +  ' -evalue 10e-5  -outfmt "6 qseqid sseqid evalue pident length bitscore ppos"'

import pandas as pd
import re
import pysnooper

#@pysnooper.snoop()
def select_blast(result1,result2,best_match=True,evalue = 10**-10, pident = 40, length = 200, bitscore = 0, ppos = 0 ,qcovs = 0):
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

if __name__ == '__main__':
    import os
    os.chdir('../../ComplementaryData/01_Sequences_analysis/blast/')
    result1 = 'Lreuteri_refseq_v01_in_Lreuteri_refseq_v02.csv'
    result2 = 'Lreuteri_refseq_v02_in_Lreuteri_refseq_v01.csv'
    result_df = select_blast(result1, result2, best_match=True,evalue = 10**-10, pident = 40, length = 0, bitscore = 0, ppos = 0)
