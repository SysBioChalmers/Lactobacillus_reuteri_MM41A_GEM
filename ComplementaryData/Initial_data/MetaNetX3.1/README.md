**MetaNetX / MNXref**

*MNXref*       |  1.0 & 1.1 |     2.0    |  3.0 & 3.1 |

*Chemicals*    |            |            |            |
BioPath        | 2010/05/03 | 2010/05/03 |            |
BiGG           | 2013/06/04 | 2015/09/02 | 2017/04/11 |
BKM/BRENDA     | 2013/04/12 |            |            |
ChEBI          |    102     |    131     |    150     |
enviPath       |            |            | 2017/04/12 |
HMDB           |    3.0     |    3.6     |    3.6     |
KEGG           | 2013/04/12 |    75.1    |    82.0    |
LipidMaps      | 2013/04/12 | 2015/06/28 | 2017/04/13 |
MetaCyc        |    17.0    |    19.1    |    20.5    |
Reactome       | 2013/05/24 |     53     |     59     |
SABIO-RK       |            |            | 2016/05/27 |
SwissLipids    |            |            | 2017/04/13 |
The SEED       | 2013/04/12 | 2013/06/19 | 2017/04/13 |
UMBBD-EAWAG    |            | 2014/06/30 |            |
UniPathway     |   2012/04  |   2015/03  |            |

*Reactions*    |            |            |            |
BioPath        | 2010/05/03 | 2010/05/03 |            |
BiGG           | 2013/06/04 | 2015/09/02 | 2017/04/11 |
BKM/BRENDA     | 2013/04/12 |            |            |
KEGG           | 2013/04/12 |    75.1    |    82.0    |
MetaCyc        |    17.0    |    19.1    |    20.5    |
Reactome       | 2013/05/24 |     53     |     59     |
Rhea           |     39     |     64     |     81     |
SABIO-RK       |            |            | 2016/05/27 |
The SEED       | 2013/04/12 | 2013/06/19 | 2017/04/13 |
UniPathway     |   2012/04  |   2015/03  |            |

*Compartments* |            |            |            |
BiGG           | 2013/06/04 | 2015/09/02 | 2017/04/11 |
CCO            |    17.0    |    19.1    |    20.5    |
GO             | 2012/06/09 | 2015/09/02 | 2017/04/11 |
The SEED       | 2013/04/12 | 2013/06/19 | 2017/04/13 |

*Proteins*
UniProt        |   2013_06  |   2015_09  |   2017_04  |


*Version 3.1*
- Incorporate new mapping fixes.
- Better identification of complex/salt chemicals.
- Merge most isotopes with their stable element.
- Force split of BiGG chemicals appearing in the same reaction.
- Add reaction original equation for most sources.
- Fix for S atom stereochemistry.
- Faster structural computations.

*Version 3.0*
- All reactions are compartmentalized now, using generic compartments.
- Integrate SABIO-RK chemicals and reactions.
- Integrate SwissLipids chemicals.
- enviPath replaces UMBBD-EAWAG.
- Do not use BioPath anymore, as our frozen data were too out-dated.
- Do not use UniPathway anymore, as most of its data are covered by Rhea and KEGG.
- Use ChEBI ontology to help merging acids/bases and tautomers.
- Add InChIKey for display, search and download.
- Compute major tautomer at pH 7.3 from SMILES or InChI now if MOL files are not
  available.


*Version 2.0*
- Sum up in "MetaNetX/MNXref - reconciliation of metabolites and biochemical
  reactions to bring together genome-scale metabolic networks", Nucleic Acids
  Research (2016) 44(D1):D523-D526. (doi: 10.1093/nar/gkv1117)
- Integrate UMBBD-EAWAG chemicals.
- Do not use BKM/BRENDA as their public exported format is more ambiguous.
- Full integration and use of MNXM01, as transport proton.
  Add a MNXR01 reaction to transform MNXM01 to and from MNXM1 (balanced proton).
- Use jChem/Marvin Beans version 6.2.0 for structural determination.
  We use the same version than Rhea to help data comparison.
- Start to use an exception file for structural determination (from Rhea)
  to avoid some jChem/Marvin Beans issues with some structures.
  We extend this exception file with ChEBI & non-ChEBI identifiers.


*Version 1.1*
- Add MNXM01, used as transport proton.
- Harmonize Unicode and non-Unicode characters to help reconciliation:
  e.g. β -> beta
- Remove fake InChI strings "unknown".
- Remove extra "\", " |" or ";" strings in some compound names to help reconciliation.
- Remove some invisible Unicode characters to help reconciliation.
- Improve some cellular compartment mapping, mainly for CCO.
- Fix MNXM114187 - ammonium hydroxide - formula.


*Version 1.0*
- Sum up in "Reconciliation of metabolites and biochemical reactions for metabolic
  networks", Brief Bioinform (2014) 15(1):123-135. (doi: 10.1093/bib/bbs058)
- Use jChem/Marvin Beans version 5.10.3 for structural determination.
