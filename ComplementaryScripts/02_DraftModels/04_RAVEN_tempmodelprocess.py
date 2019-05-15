#!/usr/bin/env python
#-*- coding: utf-8 -*-
# Created by lhao at 2019-05-07


import os
import cobra
import re
from cobra import Model, Reaction, Metabolite
from Bio.Blast import NCBIXML
import My_def



Lreu_ra_tp = cobra.io.read_sbml_model('Lreu_ra_te.xml')
'''

for i in Lreu_ra_tp.metabolites:
    if '_i' in i.id:
        for ii in i.reactions:
            ii.reaction = ii.reaction.replace(i.id,i.id.split('_i')[0])

'''


My_def.io_outtxt(Lreu_ra_tp,"Lreu_ra_te.txt",True)

