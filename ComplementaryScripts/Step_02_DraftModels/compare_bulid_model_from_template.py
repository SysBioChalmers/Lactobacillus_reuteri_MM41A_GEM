#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 2019-05-24

"""compare_bulid_model_from_template.py
:description : script
:param : 
:returns: 
:rtype: 
"""
import sys
import os



os.chdir('../../ComplementaryData/Step_02_DraftModels/Template/')
template_model = cobra.io.load_json_model('template_models/iBT721_standlized.json')

# do blast
do_blast('../Lreuteri_biogaia_v03.faa',
          'template_seqs/iBT721.faa')

BBHs = extract_BBHs(1e-20,True)

candidate_rxns = get_all_rxns_in_BBH(template_model, BBHs)

model_from_template = build_model_from_template(candidate_rxns,BBHs,True)

