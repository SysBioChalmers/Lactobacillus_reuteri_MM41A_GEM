#!/usr/bin/env python
#-*- coding: utf-8 -*-
# Created by lhao at 2019-04-05

def io_outtxt(model,outfile,sort=False):
    # function, cobrapy model write in to a txt file
    with open(outfile, 'w') as outf:
        reaidlist = [i.id for i in model.reactions ]
        if sort == True:
            reaidlist.sort()
        for reaid in reaidlist:
            rea = model.reactions.get_by_id(reaid)
            outline_list = [reaid, rea.reaction, str(rea.objective_coefficient),str(rea.lower_bound), str(rea.upper_bound),rea.gene_reaction_rule]
            outf.write("\t".join(outline_list) + '\n')
if __name__ == '__main__':
    import cobra
    import os
    os.chdir('/Users/lhao/Box Sync/Project_Lreuteri/Step_02_DraftModels/Template/')
    iNF517 = cobra.io.read_sbml_model('iNF517.xml')
    outfile = "iNF517_3.txt"
    io_outtxt(iNF517, outfile, True)
    print("Done")

