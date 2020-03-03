#!/usr/bin/env python

import urllib
from difflib import SequenceMatcher
import numpy as np
from numpy import genfromtxt
import re
import scipy.stats as st
import statsmodels.sandbox.stats.multicomp as mc

# Metabolite_Parser

def downloadParser():
    
    """
    Download compound information (cpd_id, cpd_name) from KEGG database.
    
    Returns:
        :return met2keggs: Output table with KEGG information about compound name and identifier
        :rtype met2keggs: file
    """

    met2keggs = urllib.urlopen("http://rest.kegg.jp/list/compound")
    met2keggs = met2keggs.read()
    met2keggs = met2keggs.split("\n")

    return(met2keggs)
    
def add2Dictionary(noPrefixMetName, newMetName, metabolite, kegg_cpd, keggMetNames, metDict):
    
    """
    Function to compare input metabolite name with kegg name.
     
    If newMetName (metabolite without prefix) == keggMetName:
        Calculate similarity between metabolite (oirginal name) and keggMetName
        
    It creates different dictionaries for each metabolite (one per match) 
        
     Arguments:
        :param noPrefixMetName: Input metabolite name without prefix (metabolite_parser.mainPrefixes)
        :type noPrefixMetName: string
        
        :param newMetName: Input metabolite name without prefix and chemical words.
        :type newMetName: string
        
        :param metabolite: Input metabolite name
        :type metabolite: string
        
        :param kegg_cpd: Kegg identifier of the metabolite
        :type kegg_cpd: string
        
        :param keggMetNames: "List" of Kegg metabolite names associated to a Kegg identifier
        :type keggMetNames: string
        
    Returns:
        :return metDict: Output dictionary with this structure: {Metabolite_Name_in_Kegg: [similarity, Kegg_compound_identifier]
        :rtype metDict: dictionary
    """
    
    keggMetNames = keggMetNames.split(";")
    for keggMetName in keggMetNames:
        if re.search(".*" + newMetName + ".*", keggMetName.strip(), re.IGNORECASE):
           similarity = calculateSimilarity(metabolite.strip(), keggMetName.strip())
           metDict[keggMetName.strip()] = [similarity, kegg_cpd]
           
    return(metDict)
    
def calculateSimilarity(metabolite, keggMetName):
    
    """
    Function used by add2dictionary() to compare two metabolite names and return the similarity between them.
    If the difference between them is only one of these mainPrefixes, a similarity of 90% is returned.
    
     Arguments:
        :param metabolite: Input metabolite name
        :type metabolite: string
        
        :param keggMetName: One of the metabolites contained in keggMetNames
        :type keggMetName: string
        
    Returns:
        :return similarity: Percentage of similarity between 2 input names.
        :rtype similarity: float
    """
    if (metabolite == '' or keggMetName == ''):
        raise TypeError('Empty string as input')
    mainPrefixes = {"", "cis-","trans-","d-","l-","(s)-","alpha-","beta-","alpha ","beta ","alpha-d-","beta-d-","alpha-l-","beta-l-", "l-beta-", "l-alpha-", "d-beta-", "d-alpha-"}
    if metabolite.lower() == keggMetName.lower():
        similarity = 1.0
    elif keggMetName.lower().replace(metabolite.lower(), "") in mainPrefixes:
        similarity = 0.9
    elif metabolite.lower().replace(keggMetName.lower(), "") in mainPrefixes:
        similarity = 0.9
    else:
        similarity = SequenceMatcher(a=metabolite.lower(), b=keggMetName.lower()).ratio()

    return similarity
    
# Kegg_downloader
    
def downloadKegg(args):
    
    """
    Download necessary information from Kegg Database for parsing.
    
    Returns:
        :return gen2kegg: kegg_gene_identifier "\t" Gene_Symbol ";" Gene_name
        :rtype gen2kegg: file
        
        :return kgen2pathway: kegg_gene_identifier "\t" pathway_identifier_for_gene
        :rtype kgen2pathway: file
        
        :return met2kegg: kegg_metabolite_identifier "\t" Metabolite_names_list_sep_;
        :rtype met2kegg: file
        
        :return kmet2pathway: kegg_metabolite_identifier "\t" pathway_identifier_for_metabolite
        :rtype kmet2pathway:similarity file
        
        :return pathways: pathway_identifier_for_gene "\t" Pathway_name "-" Specified_organism
        :rtype pathways: file
        
        :return pathways_common: pathway_identifier_for_metabolite "\t" Pathway_name
        :rtype pathways_common: file
    """

    urllib.urlretrieve("http://rest.kegg.jp/list/" + args.specie, args.gen2keggs)
    urllib.urlretrieve("http://rest.kegg.jp/link/pathway/" + args.specie, args.kgen2pathways)
    urllib.urlretrieve("http://rest.kegg.jp/list/compound", args.met2keggs)
    urllib.urlretrieve("http://rest.kegg.jp/link/pathway/compound", args.kmet2pathways)
    urllib.urlretrieve("http://rest.kegg.jp/list/pathway/" + args.specie, args.pathways)
    urllib.urlretrieve("http://rest.kegg.jp/list/pathway", args.pathways_commons)

    return(args)
    
def createGeneList(args):
    
    """
    Function that takes the gene expression matrix and extracts the column with Gene_Symbols
    
    Arguments:
        :param geneDataset: Gene Expression Dataset File
        :type geneDataset: file
        
        :param genId: Name of the column with Gene_symbol identifier
        :type genId: string
        
    Returns:
        :return geneList: "List" that contains only the gene symbols
        :rtype geneList: numpy_table
    """
    
    with open(args.geneDataset, "r") as data:
        header = data.readline()
	
    header = header.strip().split('\t')
    ncol = header.index(args.genId)

    geneList = genfromtxt(args.geneDataset, delimiter='\t', usecols=ncol, dtype=None)
    geneList = np.delete(geneList, 0)
    
    return(geneList)

def useMetaboliteParsed(args):
    
    """
    This function is used when the input metabolite data set has been parsed with metabolite_parser tool.
    
    Arguments:
        :param metDataset: Metabolomic Dataset File
        :type metDataset: file
        
    Returns:
        :return metList: list with selected Kegg names column of the metabolomic dataset
        :rtype metList: list
    """

    metList = []
    
    for line in open(args.metDataset, "r"):
        selected = line.split("\t")[5].strip()
        if selected == "Yes":
            metList.append(line.split("\t")[1])
    
    return(metList)
    
def useMetaboliteNotParsed(args):
    
    """
    This function is used when the input metabolite data set has NOT been parsed with metabolite_parser tool.
    
    Arguments:
        :param metDataset: Metabolomic Dataset File
        :type metDataset: file
        
        :param metId: Name of the column with metabolite names
        :type metId: string
        
    Returns:
        :return metList: "List" that contains only the metabolite names
        :rtype metList: numpy_table
    """

    # If metabolite dataset is not parsed, create a dataset with metabolite name column
    
    with open(args.metDataset, "r") as data:
        header = data.readline()
	
    header = header.strip().split('\t')
    ncol = header.index(args.metId)

    metList = genfromtxt(args.MET_dataset, delimiter='\t', usecols=ncol, dtype=None)
    metList = np.delete(metList, 0)
    
    return(metList)

def findWholeWord(word):
    
    """
    Function to find one word inside another word.
    
    Arguments:
        :param word: word to check if is contained inside another word
        :type word: string
        
    Returns:
        :return is_contained: Returns "yes" if it's contained, "no" if not
        :rtype is_contained: boolean
    """
    
    return re.compile(r'\b(?![-])({0})(?![-])\b'.format(word), flags=re.IGNORECASE).search
    
# Pathway Enrichment
    
def fisherExactTest(args, path_feat):
    
    """    
    This function performs a complete Pathway enrichment analysis:
    1) Calcule all the values of the contingency table
    2) Perform Fisher Exact test
    3) FDR correction of p-values
    
    *********************
    * Fisher exact test *
    *********************
    
             | Molecules associated to pathway | Molecules not associated to pathway | Total
    -----------------------------------------------------------------------------------------
    DEG/M    |              a                  |                 b                   |  a+b
    -----------------------------------------------------------------------------------------
    No DEG/M |              c                  |                 d                   |  c+d
    -----------------------------------------------------------------------------------------
    Total    |             a+c                 |                b+d                  |   N
    -----------------------------------------------------------------------------------------
    
    N = Total of molecules
    a+b = Total 1 flags
    c+d = Total 0 flags
    
    Arguments:       
        :params deaGeneDataset deaMetDataset: Tables with Differential Expression Analysis information for gene expression
            and metabolomics, respectively
        :types deaGeneDataset met_dataset: files
        
        :params gene_id_col, met_id_col, gene_flag_col, met_flag_col: Column names of unique identifiers and desired
            flag column of gene expression and metabolomics datasets, respectively
        :type gene_id_col, met_id_col, gene_flag_col, met_flag_col: strings
        
        :params alpha, method: alpha-value and method desired for the FDR correction
        :type alpha, method: strings
    
    Returns:
        :return output: Table with this structure: Pathway Name Odds_Ratio	P_value	FDR_Correction	Flag_#
        :rtype output: file    
    """
    
    # N)    
    N_gene = sum(1 for line in open(args.deaGeneDataset)) - 1
    N_met = sum(1 for line in open(args.deaMetDataset)) - 1
    N = N_gene + N_met
    
    # a+b, c+d)
    metDeCounter = 0
    metNonDeCounter = 0
    with open(args.deaGeneDataset, "r") as geneDataset:
        header_gene = geneDataset.readline().split("\t")
        indice_gene_flag = header_gene.index(args.gene_flag_col)
        for line in geneDataset:
            flag = int(line.split("\t")[indice_gene_flag])
            if flag == 1:
                metDeCounter += 1
            elif flag == 0:
                metNonDeCounter += 1              
               
    geneDeCounter = 0
    geneNonDeCounter = 0
    with open(args.deaMetDataset, "r") as metDataset:
        header_met = metDataset.readline().split("\t")
        indice_met_flag = header_met.index(args.met_flag_col)
        for line in metDataset:
            flag = int(line.split("\t")[indice_met_flag])
            if flag == 1:
                geneDeCounter += 1
            elif flag == 0:
                geneNonDeCounter += 1              
    
    a_b = metDeCounter + geneDeCounter
    c_d = metNonDeCounter + geneNonDeCounter
    
    # a)
    geneDataset = open(args.deaGeneDataset, "r")
    metDataset = open(args.deaMetDataset, "r")
    tmpList_v1 = []
    tmpList_v2 = []
    pvalues = []
    indice_gene_id = header_gene.index(args.gene_id_col)
    indice_met_id = header_met.index(args.met_id_col)
    
    for pathway in path_feat.keys():
        a_gen = 0
        a_met = 0
        a_c = 0
        
        for value in path_feat[pathway]:
            geneDataset.seek(0)
            geneDataset.readline()
            for line in geneDataset:
                gene = line.split("\t")[indice_gene_id].replace('"', '')
                if value == gene:
                    a_gen += int(line.split("\t")[indice_gene_flag])
                    a_c += 1
                    break
            metDataset.seek(0)
            metDataset.readline()
            for line2 in metDataset:
                metabolite = line2.split("\t")[indice_met_id].replace('"', '')
                if value == metabolite:
                    a_met += int(line2.split("\t")[indice_met_flag])
                    a_c += 1
                    break
                
        a = a_gen + a_met
        
        # b), c), d)
        b = a_b - (a_gen + a_met)
        c = a_c - (a_gen + a_met)
        d = c_d - (a_c - (a_gen + a_met))
        
        oddsratio, pvalue = st.fisher_exact([[a, b], [c, d]])
        pvalues.append(pvalue)
        
        tmpList_v1.append(pathway + "\t" + str(oddsratio) + "\t" + str(pvalue) + "\t")
    
    # FDR correction
    pvals_corrected = mc.multipletests(pvalues, alpha = float(args.alpha), method = args.method, is_sorted = False, returnsorted = False)[1]
    
    # Write results, Sort by P value
    i = 0
    for line in tmpList_v1:
        tmpList_v2.append(line + str(pvals_corrected[i]))
        i+=1
     
    with open(args.output, "w") as output:
        output.write("Pathway Name\tOdds Ratio\tP value\tFDR Correction\tFlag_0.05\n")
        for line in sorted(tmpList_v2, key=lambda line: float(line.split("\t")[2])):
            if float(line.split("\t")[3]) <= 0.05:
                output.write(line + "\t" + str(1) + "\n")
            else:
                output.write(line + "\t" + str(0) + "\n")
    
    geneDataset.close()
    metDataset.close()
    
    return(args)