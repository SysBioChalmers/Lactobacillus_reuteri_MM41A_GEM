"""
DNA
===

This module generates BOFsc for the 4 bases of DNA (dATP, dTTP, dCTP and dGTP)

"""
import gzip
import warnings

BASES = ['A', 'T', 'C', 'G']


# Methods
def _import_genome(fasta):
    from Bio import SeqIO
    try:
        # Import as a single handle genome
        extension = fasta.split('.')[-1]
        if extension == 'gz':
            with gzip.open(fasta, "rt") as handle:
                genome = list(SeqIO.parse(handle, 'fasta'))
        else:
            genome = list(SeqIO.parse(fasta, 'fasta'))
        if len(genome) > 1:
            warnings.warn(
                '%s handles in the genome file.This may indicate that your genome is not completely assembled. \nBOFdat will parse the contigs but the stoichiometric coefficients may not be accurate.' % (
                len(genome),))
    except:
        raise ImportError('The file provided cannot be imported.')

    return genome


# def _import_model(path_to_model):
#     import cobra
#     extension = path_to_model.split('.')[-1]
#     if extension == 'json':
#         return cobra.io.load_json_model(path_to_model)
#     elif extension == 'xml':
#         return cobra.io.read_sbml_model(path_to_model)
#     else:
#         raise Exception('Model format not compatible, provide xml or json')

def _get_number_of_bases(genome):
    # Counts the number of each letter in the genome
    base_genome = {'A': 0, 'T': 0, 'C': 0, 'G': 0}
    for record in genome:
        for element in record:
            value = base_genome.get(element.upper())
            if value == None:
                continue
            else:
                new_value = value + 1
                base_genome[element] = new_value

    return base_genome


def _get_ratio(base_genome):
    # Get the ratios for each letter in the genome
    ratio_genome = {'A': 0, 'T': 0, 'C': 0, 'G': 0}

    # DNA is double strand so the number of A = number of T and C = G
    at_number = base_genome.get('A') + base_genome.get('T')
    base_genome['A'] = at_number
    base_genome['T'] = at_number
    gc_number = base_genome.get('C') + base_genome.get('G')
    base_genome['G'] = gc_number
    base_genome['C'] = gc_number

    for letter in BASES:
        number_of_base = float(base_genome.get(letter))
        total = (sum(base_genome.values()))
        ratio = number_of_base / total
        ratio_genome[letter] = ratio

    return ratio_genome


def _convert_to_coefficient(ratio_genome, DNA_RATIO):
    # Transform the ratios into mmol/gDW
    DNA_WEIGHT = 1 * DNA_RATIO

    DIPHOSPHATE_WEIGHT = 174.951262

    formula_weight = {'A': 487.149863, 'T': 478.136503,
                      'C': 463.125163, 'G': 503.149263}
    base_to_bigg = {'A': 'dATP', 'T': 'dTTP',
                    'C': 'dCTP', 'G': 'dGTP'}

    # iML1515.metabolites.datp_c.formula_weight = 487.149863
    # iML1515.metabolites.dctp_c.formula_weight = 463.12516300000004
    # iML1515.metabolites.dgtp_c.formula_weight = 503.149263
    # iML1515.metabolites.dttp_c.formula_weight = 478.136503

    # base_to_bigg = {'A': model.metabolites.datp_c, 'T': model.metabolites.dttp_c,
    #                 'C': model.metabolites.dctp_c, 'G': model.metabolites.dgtp_c}

    coefficients, metabolites = [], []

    # Calculate the biomass coefficient for each metabolite
    # ratio * DNA_WEIGHT / mol_weight * 1000
    for letter in BASES:
        ratio = ratio_genome.get(letter)
        total_weight = ratio * DNA_WEIGHT
        mol_weight = formula_weight.get(letter) - DIPHOSPHATE_WEIGHT
        mmols_per_gDW = (total_weight / mol_weight) * 1000

        coefficients.append(mmols_per_gDW)
        metabolites.append(base_to_bigg.get(letter))

    DNA_coefficients = dict(zip(metabolites, [-i for i in coefficients]))
    return DNA_coefficients


def generate_coefficients(path_to_fasta, DNA_WEIGHT_FRACTION=1):
    """
    Generates a dictionary of metabolite:coefficients for the 4 DNA bases from the organism's
    DNA fasta file and the weight percentage of DNA in the cell.

    :param path_to_fasta: a path to the DNA fasta file of the organism, format should be compatible with BioPython SeqIO

    :param path_to_model: a path to the model, format supported are json and xml

    :param DNA_RATIO: the ratio of DNA in the entire cell

    :return: a dictionary of metabolites and coefficients
    """

    if DNA_WEIGHT_FRACTION > 1.:
        raise Exception('WEIGHT FRACTION should be a number between 0 and 1')
    # Operations
    genome = _import_genome(path_to_fasta)
    base_in_genome = _get_number_of_bases(genome)
    ratio_in_genome = _get_ratio(base_in_genome)

    biomass_coefficients = _convert_to_coefficient(ratio_in_genome,
                                                   DNA_WEIGHT_FRACTION)
    # Add Pyrophosphate synthesis as the sum of the coefficients
    ppi_coeff = sum(biomass_coefficients.values())
    ppi_dict = {'ppi_c': -ppi_coeff}
    biomass_coefficients.update(ppi_dict)

    return biomass_coefficients


'''
The option to update the coefficients of the metabolites in the biomass objective function is left to the user
'''
# def update_biomass_coefficients(dict_of_coefficients,model):
#     """
#     Updates the biomass coefficients given the input metabolite:coefficient dictionary.
#
#     :param dict_of_coefficients: dictionary of metabolites and coefficients
#
#     :param model: model to update
#
#     :return: The biomass objective function is updated.
#     """
#
#     from BOFdat import update
#     update.update_biomass(dict_of_coefficients,model)

if __name__ == '__main__':
    import os
    os.chdir('../../ComplementaryData/')
    # path_to_fasta = 'sequences_processing/genomic_sequences/[Bacteroides]_pectinophilus_ATCC_43243_genomic.fna.gz'
    path_to_fasta = 'Step_01_Sequences_analysis/Lreuteri_biogaia_v03/Lreuteri_biogaia_v03.fas'

    genome = _import_genome(path_to_fasta)

    base_in_genome = _get_number_of_bases(genome)
    # {'A': 2284124, 'T': 2284124, 'C': 2357528, 'G': 2357528}

    ratio_in_genome = _get_ratio(base_in_genome)

    # {'A': 0.12302322535166359,
    # 'T': 0.12302322535166359,
    # 'C': 0.1269767746483364,
    # 'G': 0.1269767746483364}

    # {'A': 0.24604645070332717,
    #  'T': 0.24604645070332717,
    #  'C': 0.2539535492966728,
    #  'G': 0.2539535492966728}

    # path_to_model = 'iML1515.json'
    # model = cobra.io.load_json_model('iML1515.json')
    biomass_coefficients = generate_coefficients(path_to_fasta, DNA_WEIGHT_FRACTION=1)
    DNA_WEIGHT_FRACTION = 0.031
    # biomass_coefficients = biomass_coefficients * DNA_WEIGHT_FRACTION
    print(biomass_coefficients)

    # DNA_coefficients = _convert_to_coefficient(model, ratio_genome, CELL_WEIGHT, DNA_RATIO)
    # path_to_model = ''
    # DNA_WEIGHT_FRACTION =
    # generate_coefficients(path_to_fasta, path_to_model , DNA_WEIGHT_FRACTION=0.031)

    # {<Metabolite dgtp_c at 0x7fdcbedefa50>: -0.023987227235417734,
    #  <Metabolite ppi_c at 0x7fdcbee83f50>: 0.10089506970390794,
    #  <Metabolite dttp_c at 0x7fdcbeebacd0>: -0.025157688898857517,
    #  <Metabolite dctp_c at 0x7fdcbeedf510>: -0.02731878216895449,
    #  <Metabolite datp_c at 0x7fdcbeef0150>: -0.024431371400678196}

    # {<Metabolite datp_c at 0x10220cc610>: -0.012215685700339098,
    #  <Metabolite dttp_c at 0x10238fd5d0>: -0.012578844449428759,
    #  <Metabolite dctp_c at 0x102392e9d0>: -0.013659391084477245,
    #  <Metabolite dgtp_c at 0x1021e69f10>: -0.011993613617708867,
    #  <Metabolite ppi_c at 0x1021e4f450>: 0.05044753485195397}

    # {'datp_c': -0.024431371400678196,
    #  'dttp_c': -0.025157688898857517,
    #  'dctp_c': -0.027318782168954493,
    #  'dgtp_c': -0.023987227235417737,
    #  'ppi_c': 0.10089506970390795}

    #  DNAS_LRE	1.37 atp_c + 0.310284595735 datp_c + 0.180906354753 dctp_c + 0.207484087701 dgtp_c + 0.301324961811 dttp_c + 1.37 h2o_c --> DNA_LRE_c + 1.37 adp_c + 1.37 h_c + 1.37 pi_c + ppi_c	0	0.0	1000.0	None	{}
