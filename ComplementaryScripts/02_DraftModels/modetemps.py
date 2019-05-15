#!/usr/bin/env python
#-*- coding: utf-8 -*-
# Created by lhao at 2019-05-02
import os
import cobra
import re
from cobra import Model, Reaction, Metabolite
from Bio.Blast import NCBIXML
import My_def


os.chdir('/Users/lhao/Box Sync/Project_Lreuteri/02_DraftModels/Template/')

iNF517 = cobra.io.read_sbml_model('../../ref/EMBL/iNF517.xml')
iNF517_bigg = cobra.io.read_sbml_model('../../ref/EMBL/iNF517_bigg.xml')

My_def.io_outtxt(iNF517,'../../ref/EMBL/iNF517.txt',True)
My_def.io_outtxt(iNF517_bigg,'../../ref/EMBL/iNF517_bigg.txt',True)

iYO844 = cobra.io.read_sbml_model('../../ref/EMBL/iYO844.xml')
iYO844_bigg = cobra.io.read_sbml_model('../../ref/EMBL/iYO844_bigg.xml')

My_def.io_outtxt(iYO844,'../../ref/EMBL/iYO844.txt',True)
My_def.io_outtxt(iYO844_bigg,'../../ref/EMBL/iYO844_bigg.txt',True)