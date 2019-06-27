#!/usr/bin/env python
#-*- coding: utf-8 -*-
# Created by lhao at 2019-04-16

import pandas as pd

def mapIDsViaMNXref(type,queryList,fromDB,toDB,dbdir = ''):
    '''
    Associate reaction/metabolite identifiers between different databases
    :param type: 'rxns' or 'mets' depending on which information is
                expected to associate
    :param queryList: (e.g. model.rxnKEGGID, model.metChEBIID)
    :param fromDB:
    :param toDB:
    :return: list
    :dbdir:
    according Hao Wang mapIDsViaMNXref.m
    example :
    targetlist, MNX_IDlist = mapIDsViaMNXref('mets',['h2o', '10fthf', 'abc', '10fthf5glu'],'bigg','seed')

    '''

    if dbdir == '':

        dbdir = '/Users/lhao/Box Sync/Projects/Project_Lreuteri/Lactobacillus_reuteri_MM41A_GEM/ComplementaryData/MetaNetX3.1/'

    if type =='rxns':
        dbfile = 'reac_xref.tsv'
    elif type == 'mets':
        dbfile = 'chem_xref.tsv'

    names = ['XREF', 'MNX_ID']

    if fromDB =='metanetx':
        querydf = pd.DataFrame({'MNX_ID': queryList, "sor": fromDB, 'order': range(len(queryList))})
    else:
        querydf = pd.DataFrame({'ids': queryList, "sor":fromDB,'order': range(len(queryList))})
    maplist = pd.read_csv(dbdir+dbfile, sep='\t',skiprows=386, names= names,usecols=[0,1])

    maplist['sor'] = maplist.XREF.apply(lambda x: x.split(':')[0])
    #maplist = maplist['sor','ids'].apply(lambda x: maplist['XREF'](x.split(':')))

    maplist = maplist[(maplist["sor"]== fromDB)|(maplist["sor"]== toDB)]
    maplist['ids'] = maplist.XREF.apply(lambda x: x.split(':')[-1])

    if fromDB =='metanetx':
        querydf_expend = pd.merge(maplist, querydf, on=['MNX_ID'], how='left').sort_values(by='order')
        MNX_IDlist = queryList

    else:
        querydf_expend = pd.merge(maplist, querydf, on=['ids','sor'], how='right').sort_values(by='order')
        querydf_expend = querydf_expend.fillna('')
        MNX_IDlist = list(querydf_expend['MNX_ID'].values)

    querydf_expend = querydf_expend.rename(columns={'ids':'fromids','sor':'fromsor'})

    if toDB =='metanetx':
        targetlist = MNX_IDlist
        return targetlist, MNX_IDlist


    target_expend = maplist[(maplist["sor"] == toDB)]
    target_expend = pd.merge(target_expend, querydf_expend, on='MNX_ID', how='right').sort_values(by='order')
    target_expend = target_expend.rename(columns={'ids':'toids','sor':'tosor'})
    target_expend = target_expend.fillna('')


    targetlist = []
    for index in range(len(queryList)):
        toids = target_expend[(target_expend['order'] == index)]['toids'].values

        if len(toids)==0:
            targele = ['']
        elif len(toids)==1:
            targele = toids[0]
        else:
            targele = list(toids)

        targetlist.append(targele)
    #targetlist = target_expend
    targetlist = list(targetlist)

    return targetlist,MNX_IDlist


    #maplistfrom = maplist[(maplist["sor"]== fromDB)]
    #maplistto = maplist[(maplist["sor"]== toDB)]

    # fromdic= (maplist["sor"]==fromDB).set_index('XREF')['MNX_ID'].to_dict()
    # todic = (maplist["sor"]==toDB).set_index('XREF')['MNX_ID'].to_dict()

if __name__ == '__main__':

    queryList =[' ',
         'MNXM237',
         'MNXR108091',
         'MNXR109012',
         'MNXR99657',
         'MNXM163072',
         'MNXR109068',
         'MNXR104316',
         'MNXR100467',
         'MNXM164047',
         'MNXR141859',
         'MNXR99471',
         'MNXR106431',
         'MNXR107867',
         'MNXR108282',
         'MNXR123144']
    #queryList = ['h2o', '10fthf', 'abc', '10fthf5glu']
    print(mapIDsViaMNXref('rxns',queryList,'metanetx','bigg'))

