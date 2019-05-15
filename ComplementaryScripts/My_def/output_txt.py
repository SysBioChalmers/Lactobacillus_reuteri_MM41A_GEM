#!/usr/bin/env python
#-*- coding: utf-8 -*-
# Created by lhao at 2019-04-23
import cobra

def output_txt(solution,outfile):
    need_fluxes =solution.fluxes[abs(solution.fluxes)>1e-10]
    with open(outfile,'w') as outf:
        for need_id in need_fluxes.index:
            rea = model.reactions.get_by_id(need_id)
            #print (need_id,need_fluxes[need_id],rea.reaction,rea.lower_bound,rea.upper_bound,rea.objective_coefficient)
            outline_list = [need_id,str(need_fluxes[need_id]),rea.reaction,str(rea.lower_bound),str(rea.upper_bound),str(rea.objective_coefficient)]
            #print (outline_list)
            outf.write("\t".join(outline_list)+'\n')
if __name__ =="__main__":
    infile = 'D:/Work/1/textbook.xml'
    outfile = 'D:/Work/1/textbook_result.txt'

    model = cobra.io.read_sbml_model('D:/Work/1/textbook.xml')
    solution = model.optimize()
    # solution = cobra.flux_analysis.pfba(model)
    output_txt(solution,outfile)
    print ('done')