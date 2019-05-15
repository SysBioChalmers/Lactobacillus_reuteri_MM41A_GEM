#!/usr/bin/env python
#-*- coding: utf-8 -*-
# Created by lhao at 2019-04-20

from matplotlib_venn import venn3, venn3_circles
from matplotlib import pyplot as plt

def venn3_samesize(subsets, set_labels=('A', 'B', 'C'), set_colors=('r', 'g', 'b'), alpha=0.4, normalize_to=1.0, ax=None, subset_label_formatter=None):

    v = venn3(subsets=(1, 1, 1, 1, 1, 1, 1),set_labels=set_labels,
              set_colors=set_colors, alpha=alpha, normalize_to=normalize_to,
              ax=ax, subset_label_formatter=subset_label_formatter)

    a = subsets[0]
    b = subsets[1]
    c = subsets[2]
    listid = ['100', '010', '110', '001', '101', '011', '111']
    testlist = [
    len(a - (b | c)),  # TODO: This is certainly not the most efficient way to compute.
    len(b - (a | c)),
    len((a & b) - c),
    len(c - (a | b)),
    len((a & c) - b),
    len((b & c) - a),
    len(a & b & c)]
    #testlist = [str(i) for i in testlist ]
    for i in range(7):
        v.get_label_by_id(listid[i]).set_text(testlist[i])

    return v
