#!/usr/bin/env python
#-*- coding: utf-8 -*-
# Created by lhao at 2019-05-15
#    find BBH(Bidirectional Best Hits) from blast results
#     balstcmd = 'blastp -db ' + db + ' -query ' + seq + ' -out ' + outfile +  ' -evalue 10e-5  -outfmt "6 qseqid sseqid evalue pident length bitscore ppos"'

import pandas as pd
import re
import pysnooper
import My_def

@pysnooper.snoop()
def select_blast(result1,result2,best_match=True,evalue = 10^-10, pident = 40, length = 200, bitscore = 0, ppos = 0 ):
    names = ['qseqid', 'sseqid', 'evalue', 'pident', 'length', 'bitscore', 'ppos']

    df1 = pd.read_csv(result1, sep = '\t', names = names)
    df2 = pd.read_csv(result2, sep = '\t', names = names)
    
    for index  in [0,1]:
        dfi = [df1,df2][index]
        dfi = dfi[(dfi["evalue"] <= evalue) & (dfi["pident"] >= pident) & (dfi["length"] >= length) & (
                dfi["bitscore"] >= bitscore) & (dfi["ppos"] >= ppos)]
        dfi['sseqid'] = dfi.sseqid.apply(lambda x: re.sub('(\|$)|(^ref\|)','',x))
        dfi = dfi.sort_values(names, ascending=[False,False,False,True,True,True,True])

        if best_match:
            dfi = dfi.drop_duplicates(subset='qseqid', keep='first')
        if index ==1:
            df1 = dfi
        else:
            df2 = dfi
    df2 = df2.rename(columns={'qseqid':'sseqid','sseqid':'qseqid'})

    result_df = pd.merge(df1, df2[['qseqid','sseqid']], on=['qseqid','sseqid'], how='inner')
    return df1,df2,result_df
if __name__ == '__main__':
    import os
    from Bio.Blast import NCBIXML

    os.chdir('../../ComplementaryData/01_Sequences_analysis/blast/')
    result1 = 'Lreuteri_refseq_v01_in_Lreuteri_refseq_v02.csv'
    result2 = 'Lreuteri_refseq_v02_in_Lreuteri_refseq_v01.csv'
    df1,df2,result_df = select_blast(result1, result2, best_match=True)
