import os
import cobra
import re
from Bio import SeqIO


def do_blast(query_fasta,target_fasta):
    # query_fasta: fasta file for the organsims to be modelled
    # target_fasta: fasata file that contains enzyme sequences in the templates model
    
    # combine two fasta files, do blast. Use default evalue. The blast results will be filetered based on evalue
    # in the BBH extraction step
    
    # build combined fasta file for blast
    fhand = open('input.fa','w')
    for rec in SeqIO.parse(query_fasta,'fasta'):
        fhand.write('>query|{0}\n{1}\n'.format(rec.id,rec.seq))
    
    for rec in SeqIO.parse(target_fasta,'fasta'):
        fhand.write('>target|{0}\n{1}\n'.format(rec.id,rec.seq))
    
    fhand.close()
    
    blast_cmd = '''mkdir blast_tmps
    mv input.fa blast_tmps
    makeblastdb -dbtype prot -in blast_tmps/input.fa  -out blast_tmps/DB
    blastp -query blast_tmps/input.fa -db blast_tmps/DB -outfmt 6  -max_target_seqs 500 -max_hsps 1  -num_threads 1 -out blast/bidirectional_blast.tab
    '''
    os.system(blast_cmd)
    add_coverage('blast_tmps/input.fa','blast/bidirectional_blast.tab')
    os.system('rm -r blast_tmps')


def add_coverage(combined_fasta,blast_results):
    # combined_fasta: combined fasta file, which contains all sequnces from both query and target proteome.
    # blast_results: tsv file produced by blastp
    # 
    # Output: add a new column: coverage at the end of the blast_results.
    # Coverage is caculated as: length of aligned query sequence/length of query sequence * 100
    # 
    
    # Load sequence lenght
    Lenghs = dict()
    for rec in SeqIO.parse(combined_fasta,'fasta'): Lenghs[rec.id] = len(rec.seq)
    
    
    # calculate coverage
    outfile = blast_results.replace('.tab','_cov.tab')
    if outfile == blast_results: outfile += '_cov.tab'
    
    fhand = open(outfile,'w')
    for line in open(blast_results):
        cont = line.split()
        aln_lengh = float(cont[7])-float(cont[6])
        cov = aln_lengh/Lenghs[cont[0]]*100
        fhand.write(line.strip()+'\t{0}\n'.format(cov))
    fhand.close()
    
    # replace the orginal blast results file with the new one.
    os.system('mv {0} {1}'.format(outfile,blast_results))


def extract_BBHs(BBH_evalue,report=False,coverage=45):
    # BBH_evalue: like 1e-10
    # report: if True, print out the nubmer of BBHs
    # coverage: the lower bound of query coverage to determine a true hit. from 0-100
    # 
    
    target_query = dict()
    query_target = dict()
    for line in open('blast/bidirectional_blast.tab'):
        cont = line.split()
        evalue = float(cont[10])
        cov = float(cont[12])
        # filter based on evalue and coverage
        if evalue > BBH_evalue: continue
        if cov < coverage: continue
        
        if cont[0].startswith('target|') and cont[1].startswith('query|'):
            target_gene = cont[0].split('|')[1]
            query_gene = cont[1].split('|')[1]
            if target_query.get(target_gene,(None,1000))[1]>evalue: target_query[target_gene] = (query_gene,evalue)
        elif cont[0].startswith('query|') and cont[1].startswith('target|'):
            target_gene = cont[1].split('|')[1]
            query_gene = cont[0].split('|')[1]
            if query_target.get(query_gene,(None,1000))[1]>evalue: query_target[query_gene] = (target_gene,evalue)
        else: None

    # 3. find BBH, the overlaped gene pairs
    BBH = list()
    for target_gene, (query_gene,evalue) in target_query.items():
        if query_target.get(query_gene,(None,1000))[0] == target_gene: BBH.append((query_gene,target_gene))
    
    if report:
        print('Number of homologs for query genes (evalue = {0}):'.format(BBH_evalue),len(query_target))
        print('Number of homologs for target genes (evalue = {0}):'.format(BBH_evalue),len(target_query))
        print('Number of BBHs (evalue = {0}):'.format(BBH_evalue),len(BBH))
    return BBH


def get_all_rxns_in_BBH(template_model, BBHs):
    # model file: templates model file
    # BBH, best bidirectional hit

    # get all reactions
    candidate_rxn_ids = list()
    candidate_rxns = list()
    gene_set = set([ i.id for i in template_model.genes])
    for query_gene,target_gene in BBHs:
        if target_gene in gene_set:     # KeyError  KeyError: 'lp_0001' (gene id not in model)
            gene = template_model.genes.get_by_id(target_gene)
            for rxn in gene.reactions:
                if rxn.id not in candidate_rxn_ids:
                    candidate_rxn_ids.append(rxn.id)
                    candidate_rxns.append(rxn)
    print('Number of candiate reactions:',len(candidate_rxn_ids))
    
    return candidate_rxns


def build_model_from_template(candidate_rxns,BBHs,report=False):
    rxns_to_add = list()
    for rxn in candidate_rxns:
        new_gr = update_gr(rxn.gene_reaction_rule,BBHs)
        new_rxn = rxn.copy()
        new_rxn.gene_reaction_rule = new_gr
        rxns_to_add.append(new_rxn)
    model_from_template = cobra.Model('model_from_template')
    model_from_template.add_reactions(rxns_to_add)
    if report: report_model_status(model_from_template)
    return model_from_template


def parse_gr(gr):
    # return a list of genes in gr rule
    genes = list(set([item for item in gr.replace('(','').replace(')','').split() if item not in ['and','or']]))
    genes.sort()
    return genes


def update_gr(gr,BBHs):
    # replace the targetli gene ids in gr rule with the ones in querymonas based on BBH. For those targetli genes
    # wihout BBH in querymonas, add a tag '_missing' at the end of the gene id
    
    target2query_map = dict()
    for query_gene,target_gene in BBHs: target2query_map[target_gene] = query_gene
    
    genes = parse_gr(gr)
    new_gr = gr
    for gene in genes:
        query_gene = target2query_map.get(gene)
        if query_gene is None: query_gene = gene+'_missing'
        new_gr = new_gr.replace(gene,query_gene)
    
    # refine new gr rule. 
     
    # case 1: if 'missing' not in gr rule, do not need to update
    if '_missing' not in new_gr: refined_gr = new_gr
    
    # case 2: only 'or' relationships are in gr rule, remove all genes with '_missing' tag
    elif ' and ' not in new_gr: 
        new_genes = [item for item in parse_gr(new_gr) if '_missing' not in item]
        refined_gr = ''
        for g in new_genes: refined_gr += g + ' or '
        refined_gr = refined_gr[:-4]
    
    # case 3: if 'or' not in gr rule, keep missing genes
    elif 'or' not in new_gr: refined_gr = new_gr
    
    # case 4: if (A and B) or (C and D_missing) or E. In such or structure, if one of the gene/gene complex is 
    #         complete, then remove all others with missing gene(s). if all parts are incompelte, only removes 
    #         those part which are totally missing.
    #elif 'and' in new_gr and 'or' in new_gr: refined_gr = update_and_or_case(new_gr)

    # Ohter cases: export to a txt file and manully edit
    else:  refined_gr = new_gr
    return refined_gr


def report_model_status(model):
    print('Number of reactions:',len(model.reactions))
    print('Number of metabolits:',len(model.metabolites))
    print('Number of compartments:',len(model.compartments))
    k = 0
    for gene in model.genes:
        if 'missing' in gene.id: k += 1
    print('Number of genes:',len(model.genes))
    print('Number of missing genes:',k)
    
    k = 0
    for rxn in model.reactions:
        if '_missing' in rxn.gene_reaction_rule: k += 1
    print('Number of reactions with missing genes:',k)


if __name__ =='__main__':
    # Modeling pipeline

    # load templates model
    #template_model = cobra.io.read_sbml_model('iJO1366.xml')
    os.chdir('../../ComplementaryData/Step_02_DraftModels/Template/')
    template_model = cobra.io.load_json_model('template_models/iBT721_standlized.json')

    # do blast
    do_blast('../Lreuteri_biogaia_v03.faa',
              'template_seqs/iBT721.faa')

    BBHs = extract_BBHs(1e-20,True)

    candidate_rxns = get_all_rxns_in_BBH(template_model, BBHs)

    model_from_template = build_model_from_template(candidate_rxns,BBHs,True)

