#!/usr/bin/env python
#-*- coding: utf-8 -*-
# Created by lhao at 2019-04-12

import cobra
import re
def gpr2log( gpr_i, othergeneset,inset = True, notinset = False, emptyset = '' ):
    torf = ''
    if  not gpr_i :
        if emptyset == '' :
            new_gpr_i = gpr_i
        else:
            new_gpr_i = emptyset
    else:
        _, geneset = cobra.core.gene.parse_gpr(gpr_i)
        new_gpr_i = gpr_i
        for gen_i in geneset:
            gen_i_re = re.compile('(^|(?<=( |\()))%s(?=( |\)|$))' %
                                  re.escape(gen_i))

            if gen_i in othergeneset:
                temp_gen_i = inset
            else:
                temp_gen_i = notinset
            new_gpr_i = gen_i_re.sub(str(temp_gen_i),new_gpr_i)


        '''for gen_i in geneset:
            if gen_i in othergeneset:
                temp_gen_i = inset
            else:
                temp_gen_i = notinset
            new_gpr_i = new_gpr_i.replace(gen_i, str(temp_gen_i))
        '''


        if (inset == True and notinset == False):
            torf = eval(new_gpr_i)
    return new_gpr_i,torf

if __name__ == '__main__':
    import cobra
    import cobra.test
    data_dir = cobra.test.data_dir
    textbook_model = cobra.test.create_test_model("textbook")
    allgenesset = [i.id for i in textbook_model.genes[0:-10]]
    gprlist = [i.gene_reaction_rule for i in textbook_model.reactions]
    for gpr_i in gprlist:
        new_gpr_i,torf = gpr2log( gpr_i ,allgenesset )
        print(gpr_i, "\t", new_gpr_i,'\t', torf)

    print("Done")

