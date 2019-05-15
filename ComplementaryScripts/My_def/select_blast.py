#!/usr/bin/env python
#-*- coding: utf-8 -*-
# Created by lhao at 2019-05-15
#    find BBH(Bidirectional Best Hits) from blast results
#     balstcmd = 'blastp -db ' + db + ' -query ' + seq + ' -out ' + outfile +  ' -evalue 10e-5  -outfmt "6 qseqid sseqid evalue pident length bitscore ppos"'

from Bio.Blast import NCBIXML
import pandas as pd

def select_blast(result1,result2,best_match=True,evalue = 10^-10, pident = 40, length = 200, bitscore = 0, ppos = 0 ):
    names = ['qseqid', 'sseqid', 'evalue', 'pident', 'length', 'bitscore', 'ppos']

    df1 = pd.read_csv(result1,names = names)
    df2 = pd.read_csv(result2, names=names)
    
    for pd_i  in [df1,df2]:
        dfi = dfi[(dfi["evalue"] <= evalue) | (dfi["pident"] >= pident) | (dfi["length"] >= length) | (
                dfi["bitscore"] >= bitscore) | (dfi["ppos"] >= ppos)]
        dfi = dfi.sort_values('evalue', ascending=False)

        if best_match:
            dfi = dfi.drop_duplicates(subset='qseqid', keep='first')


if __name__ == '__main__':
    import os
    from Bio.Blast import NCBIXML

    os.chdir('../../ComplementaryData/01_Sequences_analysis/')
    result1 = 'Lreuteri_refseq_v01_in_Lreuteri_refseq_v02.csv'
    result2 = 'Lreuteri_refseq_v01_in_Lreuteri_refseq_v02.csv'
    select_blast(result1, result2, best=True, bitscore=100, ppos=45)