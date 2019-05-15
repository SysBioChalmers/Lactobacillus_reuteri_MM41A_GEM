#!/usr/bin/env python
#-*- coding: utf-8 -*-
# Created by lhao at 2019-04-12
__all__ = ['eva_gprinset']

import cobra
def eva_gprinset( reaction , othergeneset, emptyset = '' ):

    torf = ''
    if not reaction.genes:
        return emptyset
    for gene_i in reaction.genes:
        gene_i.functional = False
        if gene_i.id in othergeneset:
            gene_i.functional = True
    return reaction.functional

if __name__ == '__main__':
    import cobra.test
    data_dir = cobra.test.data_dir
    textbook_model = cobra.test.create_test_model("textbook")
    allgenesset = [i.id for i in textbook_model.genes[0:-10]]
    new_gprlist = []
    for reaction in textbook_model.reactions:
        oldfunctioal = reaction.functional
        torf = eva_gprinset( reaction ,allgenesset )
        print(reaction.gene_reaction_rule, '\t', oldfunctioal,'\t', reaction.functional)

    print("Done")