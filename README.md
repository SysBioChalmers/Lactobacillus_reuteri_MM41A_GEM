# Lactobacillus_reuteri_MM41A_GEM
GEM for Lactobacillus reuteri ATCC PTA 6475/MM4-1A

Requirements,
- [Prokka](http://github.com) 
- [Blast+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
- [dimon](https://github.com/bbuchfink/diamond)
- [Biopython](https://biopython.org/)
- [cobrapy](https://opencobra.github.io/cobrapy/)

Notice, <br />
Please add `./ComplementaryScripts/My_def` to sys.path to ues functions in My_def   <br />
    - `import sys`    <br />
    - `sys.path.extend(['./Projects/Project_Lreuteri/Lactobacillus_reuteri_MM41A_GEM/ComplementaryScripts/My_def'])` 

Structure,

- ComplementaryData/  <br />
    Initial data/  <br />
    - Sequences of *L.reuteri* ,
        - Lreuteri_refseq_v01, Genbank , [FTP ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_refseq/ASSEMBLY_BACTERIA/Lactobacillus_reuteri/GCF_000159475/](ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_refseq/ASSEMBLY_BACTERIA/Lactobacillus_reuteri/GCF_000159475/)
        - Lreuteri_refseq_v02, Genbank , ACCESSION   [NZ_ACGX02000007,FTP ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/159/475/GCF_000159475.2_ASM15947v2/](ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/159/475/GCF_000159475.2_ASM15947v2/)
        - Lreuteri_biogaia_v03, Biogaia company, raw genome Lreuteri_biogaia_v03_Complete_sequence.fas annotated by [Prokka](http://github.com) 
    - Template models informations,
        - [iNF517](http://bigg.ucsd.edu/models/iNF517 ), Lactococcus lactisMG1363, BiGG
            - seq, GenBank, GCF_000009425.1
            - 2006, 
            - Journal of Biological Chemistry	
            - Analysis of growth of Lactobacillus plantarum WCFS1 on a complex medium using a genome-scale metabolic model
            
        - [iBT721](https://www.ebi.ac.uk/biomodels/MODEL1507180045 ), Lactobacillus plantarumWCFS1, BioModels, 
            - seq, GenBank, GCF_000203855.3
            - 2013, 
            - Applied Microbiology and Biotechnology,	
            - Genome-scale metabolic model for Lactococcus lactisMG1363 and its application to the analysis of flavor formation
        - [iML1515](http://bigg.ucsd.edu/models/iML1515 ), Escherichia coli str. K-12 substr. MG1655, , BiGG
            - seq, GenBank, [NC_000913.3](http://bigg.ucsd.edu/genomes/ncbi_accession:NC_000913.3 )
            - 2017
            - Nature Biotechnology
            - [iML1515, a knowledgebase that computes Escherichia coli traits](10.1038/nbt.3956 )
            
        - Potential template models,
          Streptococcus thermophilus LMG18311	iMP429	biomodels,
          Bacillus subtilis subsp. subtilis str. 168	iYO844	bigg,

    - Database,
        - [BiGG](http://bigg.ucsd.edu/data_access)
        - [MetaNetX](https://www.metanetx.org/mnxdoc/mnxref.html)
   interinal_data
   
 - ComplementaryScripts/

     - 01_Sequences_analysis
     - 02_DraftModels
     - 03_Compare_Refine
     - 04_Simulation
 









