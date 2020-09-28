from biopandas.pdb import PandasPdb
import math
import pandas as pd
import numpy as np
from scipy.spatial.distance import pdist
import itertools
import argparse
import os
from io import StringIO
import urllib.parse
import urllib.request
from Bio import SeqIO
from Bio.Alphabet import generic_protein
from Bio import Align

protein_letters_3to1 = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K', 'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

### IO FUNCTIONS

# Set up command-line arguments
def get_args():
    parser = argparse.ArgumentParser(
        description='Run phosphosite distance analysis')
    parser.add_argument('phosphosites', type=str)
    parser.add_argument('--path', type=str)
    parser.add_argument('--rcsb', type=str)
    return parser.parse_args()

# Uses four-letter pdb code to get .pdb file from rcsb.org. Returns an error if invalid code.
def get_rcsb(pdb_id):
    print("RCSB...", pdb_id)
    try:
        ppdb = PandasPdb().fetch_pdb(pdb_id)
        print("Fetched", pdb_id)
        return ppdb
    except:
        print("PDB file not found on RCSB", pdb_id)
        raise SystemExit

# Gets the PDB file from local user
def get_file_input(path):
    print("Local file input...")
    try:
        ppdb = PandasPdb().read_pdb(path)
        print("Read")
        return ppdb
    except:
        print("No PDB file found at " + path)
        raise SystemExit

# Determines where to get pdb file from RCSB or local drive
def get_ppdb(args):
    if args.path:
        return get_file_input(args.path)
    elif args.rcsb:
        return get_rcsb(args.rcsb)
    else:
        print("No input file or PDB code")
        raise SystemExit



# Gets phosphosite list from local drive

def get_phosphosites(path):
    if not os.path.exists(path):
        print("No phosphosites found")
        return pd.Series()
    try:
        print(path)
        print( pd.read_csv(path, squeeze=True))
        return pd.read_csv(path, squeeze=True)
        '''
        df_split = split_by_chain(pd.read_csv(path, names=['chain_id', 'residue_number'], header=0))
        phosphosites = {}
        for i in df_split.groups.keys():
            phosphosites[i] = np.unique(df_split.get_group(i).loc[:,'residue_number'])
        return phosphosites
        '''
    except:
        print("No CSV file found for phosphosites")
        raise SystemExit

# Runs the input process
def get_input():
    args = get_args()
    ppdb = get_ppdb(args)
    phosphosites = get_phosphosites(args.phosphosites)
    return ppdb, phosphosites

# Writes the final dataframe to an output file
def write_df(df, file_name, dir = "./phosphosite_analysis_results"):
    path = dir + "/" + file_name
    if not os.path.exists(dir):
        os.mkdir(dir)
    if os.path.exists(path):
        os.remove(path)
    df.to_csv(path, index=False)

def write_phos(phosphosites, uniprot_id, dir):
    path = dir + "/" + uniprot_id + "_phos.csv"
    if not os.path.exists(dir):
        os.mkdir(dir)
    if os.path.exists(path):
        os.remove(path)
    phosphosites.to_csv(path, header=False, index=False)

### DISTANCE FUNCTIONS
# Filters the atoms to keep only the atoms that are on the appropriate residues and are the correct oxygens
def get_phospho_coords_chain(residue_nums, df):
    residue_criteria = df['residue_number'].isin(residue_nums)
    serine_criteria = (df['residue_name'] == 'SER') & (df['atom_name'] == 'OG')
    threonine_criteria = (df['residue_name'] == 'THR') & (
        df['atom_name'] == 'OG1')
    tyrosine_criteria = (df['residue_name'] == 'TYR') & (
        df['atom_name'] == 'OH')
    alt_loc_criteria = (df['alt_loc'] == '') | (df['alt_loc'] == 'A')
    return df[residue_criteria & (serine_criteria | threonine_criteria | tyrosine_criteria) & alt_loc_criteria]

# Applies get_phospho_coords_chain() on all chains of the protein
def get_phospho_coords_protein(df_split, phosphosites):
    df = pd.DataFrame()
    for chain in df_split.groups.keys():
        if chain in phosphosites:
            df = pd.concat([df, get_phospho_coords_chain(
                phosphosites[chain], df_split.get_group(chain))])
    return df

# Gives pairwise combinations of all residues
def get_residue_pairs_euclidean(df):
    return np.asarray(list(itertools.combinations(df['residue_number'].tolist(), 2))), np.asarray(list(itertools.combinations(df['chain_id'].tolist(), 2)))

# Gives pairwise euclidean distances of all input residues
def get_euclidean_distances(df):
    coords = df.iloc[:, 11:14]
    return pdist(coords.to_numpy())

# Gives pairwise combination of all possible linear distances on a chain
def get_residue_pairs_linear_chain(phosphosites, chain):
    y = np.hsplit(np.asarray(list(itertools.combinations(phosphosites, 2))), 2)
    df = pd.DataFrame.from_dict(
        {"Residue A": np.squeeze(y[0]), "Residue B": np.squeeze(y[1])})
    df.insert(0, 'Chain A', chain)
    df.insert(2, 'Chain B', chain)
    return df

# Applies get_residue_pairs_linear_chain() to the whole protein
def get_residue_pairs_linear_protein(df_split, phosphosites):
    df = pd.DataFrame(columns=['Chain A', 'Residue A', 'Chain B', 'Residue B'])
    for i in df_split.groups.keys():
        if i in phosphosites:
            df = pd.concat([df, get_residue_pairs_linear_chain(phosphosites[i], i)])
    return df

### DATA EXTRACTION METHODS

# Gets dbref entry from PDB
def get_dbref(ppdb):
    others = ppdb.df['OTHERS']
    dbref = others[others['record_name'] == 'DBREF']
    return dbref['entry']

# Gets seqres entry from PDB
def get_seqres(ppdb):
    others = ppdb.df['OTHERS']
    seqres = others[others['record_name'] == 'SEQRES']
    return seqres['entry']

# Uses dbref to create a dictionary matching each chain to a gene (represented by UniProt ID)
def get_chain_gene(ppdb):
    entry = get_dbref(ppdb)
    chain_gene_dict = {}

    for i in entry:
        lst = i.split()
        if lst[4] != 'UNP':
            print("Not UniProt. Error.")
        chain_gene_dict[lst[1]] = [
            lst[5], np.arange(int(lst[2]), int(lst[3]) + 1)]

    return pd.DataFrame.from_dict(chain_gene_dict, orient='index', columns=['gene ID', 'amino acids'])

# Get a dictionary of
def get_chain_gene_dict(ppdb):
    chain_gene = get_chain_gene(ppdb)
    return chain_gene.iloc[:,:1].squeeze().to_dict()

# Splits a protein by chain
def split_by_chain(df):
    return df.groupby(['chain_id'])

# Gets header
def get_header(ppdb):
    others = ppdb.df['OTHERS']
    header = others[others['record_name'] == 'HEADER']
    return header['entry']

# Gets pdb code from header
def get_pdb_code(ppdb):
    return get_header(ppdb)[0].split()[-1]
    
### SGD METHODS

# Returns a list of amino acids by parsing the SEQRES region
def get_seqres_chain(seqres):
    amino_acids = []
    for i in seqres:
        for j in i.split()[3:]:
            amino_acids.append(protein_letters_3to1[j])
    return ''.join(amino_acids)

# Gets the SGD ID from a uniprot ID through the UniProt database
def get_sgd_id(uniprot_id):
    url = 'https://www.uniprot.org/uploadlists/' 
    
    params = {
    'from': 'ACC+ID',
    'to': 'SGD_ID',
    'format': 'tab',
    'query': uniprot_id
    }

    data = urllib.parse.urlencode(params)
    data = data.encode('utf-8')
    req = urllib.request.Request(url, data)
    with urllib.request.urlopen(req) as f:
        response = f.read()
    return (response.decode('utf-8').split()[2:][1])

# Uses aligner function to check the percentage the SEQRES sequence matches the SGD sequence
def check_sgd_chain(seqres, sgd):
    score = Align.PairwiseAligner().score(seqres, sgd)
    return score/float(len(sgd)), score/float(len(seqres))

# Runs check_sgd_chain() on the whole protein
def check_sgd_protein(ppdb, chain_gene, fasta_dict, sgd_mapping):
    
    seqres = get_seqres(ppdb)
    match_val_dict = {}

    for i in chain_gene.index.values:
        uniprot = chain_gene.loc[i, 'gene ID']
        if uniprot in sgd_mapping:
            print(sgd_mapping[uniprot])
            seq = get_seqres_chain([elem for elem in get_seqres(ppdb) if elem.split()[1] == i])
            match_val_sgd, match_val_seqres = check_sgd_chain(seq, fasta_dict[sgd_mapping[uniprot]])
            match_val_dict[i] = (match_val_sgd, match_val_seqres)
    return match_val_dict

# Generates a dictionary from the SGD data mapping each SGD ID to a protein sequence
def get_fasta_dict():
    fasta_dict = {}
    for seq_record in SeqIO.parse("../orf_trans.fasta", "fasta", generic_protein):
        fasta_dict[seq_record.description.split()[2][6:-1]] = seq_record.seq[:-1]
    return fasta_dict
    
### ERROR-CHECKING METHODS

# Determines which phosphosites are missing by checking if they do not appear in the 'ATOM' section
def phosphosites_missing(df_split, phosphosites):
    missing = {}
    for i in df_split.groups.keys():
        if i in phosphosites:
            missing[i] = [elem for elem in phosphosites[i]
                          if elem not in df_split.get_group(i)['residue_number'].values]
    return missing

# Looks at the sequence of a chain and tells which amino acids are missing or extra compared to reference sequence
def validate_chain_sequence(df, amino_acids):

    pdb_residue_nums = np.unique(df['residue_number'].to_numpy())

    extra = np.setdiff1d(pdb_residue_nums, amino_acids)

    print("Extra: ", extra)
    missing = np.setdiff1d(amino_acids, pdb_residue_nums)
    print("Missing: ", missing)
    return extra, missing

# Finds the discrepancy in phosphosites between a chain and a reference sequence 
# i.e. phosphosites that should be there according to DBREF but are not there

def validate_protein_sequence(df_split, chain_gene, missing_phosphosites):
    amino_acids = chain_gene.loc[:, 'amino acids']
    missing = {}
    discrepancy = {}
    for i in df_split.groups.keys():
        print("\n\nChain: " + i)
        missing[i] = validate_chain_sequence(
            df_split.get_group(i), amino_acids.loc[i])[1]
        if i in missing_phosphosites:
            discrepancy[i] = np.setdiff1d(missing_phosphosites[i], missing[i])
    return discrepancy

# Analyzes the discrepancy to see if there is an error
def analyze_discrepancy(discrepancy):
    error = False
    print("\n\n")
    for i in discrepancy:
        if discrepancy[i].size > 0:
            error = True
            print("Discrepancy in Chain", i, ":", discrepancy[i])
    if error:
        return "Error!"
    else:
        return "No error!"

def out_of_range(phosphosites, chain_gene):
    dic = {}
    print(phosphosites)
    for chain in phosphosites:
        phos = phosphosites[chain]
        supposed_range = chain_gene.loc[chain, 'amino acids']
        dic[chain] = any(item in phos for item in supposed_range)
    return dic

def add_surface_residues(pdb_id, df):
    surface_residues = pd.read_csv('../surface_residues/' + pdb_id + '_surfaceresidues.csv')
    print(surface_residues)

    df_A = pd.merge(df, surface_residues, left_on = ['Chain A', 'Residue A'], right_on = ['Chain', 'Residue'], how='inner')
    df_B = pd.merge(df, surface_residues, left_on = ['Chain B', 'Residue B'], right_on = ['Chain', 'Residue'], how='inner')

    tuple_list_A = list(df_A.loc[:,['Chain A', 'Residue A']].to_records(index=False))
    tuple_list_B = list(df_B.loc[:,['Chain B', 'Residue B']].to_records(index=False))
    
    tuple_list_A = [elem.tolist() for elem in tuple_list_A]
    tuple_list_B = [elem.tolist() for elem in tuple_list_B]

    df['Surface A'] = df.apply(lambda elem: (elem['Chain A'], elem['Residue A']) in tuple_list_A, axis=1)
    df['Surface B'] = df.apply(lambda elem: (elem['Chain B'], elem['Residue B']) in tuple_list_B, axis=1)

    return df

def add_interface_residues(df, interface_dict):
    pass
    
### MAIN METHOD
def analyze_phosphosite_distances(pdb_id, ppdb, phosphosites, uniprot_common, uniprot_systematic, sgd_mapping, fasta_dict): #interface_residues):
    
    # Analyzes the file to map each chain to a gene
    chain_gene = get_chain_gene(ppdb)
    chain_gene_dict = chain_gene.iloc[:,:1].squeeze().to_dict()

    print('here')

    #print(sgd_mapping)
    print(check_sgd_protein(ppdb, chain_gene, fasta_dict, sgd_mapping))

    print('hereeee')
    print(phosphosites)
    print(out_of_range(phosphosites, chain_gene))
    # Sorts phosphosites in each chain
    for elem in phosphosites:
        phosphosites[elem].sort()

    print('here2')

    # Finds phosphosites in the protein as well as the missing phosphosites
    df = ppdb.df['ATOM']
    df_filtered = get_phospho_coords_protein(split_by_chain(df), phosphosites)
    missing_phosphosites = phosphosites_missing(
        split_by_chain(df_filtered), phosphosites)

    print('here3')

    # Gets all linear pairs of residues
    linear_pairs = get_residue_pairs_linear_protein(split_by_chain(df), phosphosites)

    # Validates the sequence by analyzing discrepancy
    # i.e. if a phosphosite is missing when it should not be
    discrepancy = validate_protein_sequence(
        split_by_chain(df), chain_gene, missing_phosphosites)
    print(analyze_discrepancy(discrepancy))

    print('here4')
    # Processes the euclidean distance pairs
    residue_pairs, chain_pairs = get_residue_pairs_euclidean(df_filtered)
    residue_pairs_split, chain_pairs_split = np.hsplit(residue_pairs, 2), np.hsplit(chain_pairs, 2)
    residue_1, chain_1 = residue_pairs_split[0], chain_pairs_split[0]
    residue_2, chain_2 = residue_pairs_split[1], chain_pairs_split[1]

    # Gets euclidean distances
    euclidean_distances = get_euclidean_distances(df_filtered)

    print('here5')

    # Creates a final dataframe with the two residues and their chains and euclidean_distances
    dictionary = {"Chain A": np.squeeze(chain_1), "Residue A": np.squeeze(residue_1), "Chain B": np.squeeze(
        chain_2), "Residue B": np.squeeze(residue_2), "3D distance": euclidean_distances}
    df_final = pd.DataFrame.from_dict(dictionary)

    # Merges the linear distance pairs such that only the pairs that don't already exist in the euclidean section are appended
    df_final = df_final.merge(linear_pairs, how='outer', left_on=['Chain A', 'Residue A', 'Chain B', 'Residue B'], right_on=[
                              'Chain A', 'Residue A', 'Chain B', 'Residue B'])

    # Finds linear distance between the residues
    df_final['Linear Distances'] = np.absolute(np.subtract(df_final['Residue A'], df_final['Residue B']))

    df_final['UniProt A'] = df_final.loc[:, 'Chain A'].apply(lambda elem: chain_gene_dict[elem])
    df_final['UniProt B'] = df_final.loc[:, 'Chain B'].apply(lambda elem: chain_gene_dict[elem])
    df_final['Common Gene A'] = df_final.loc[:, 'UniProt A'].apply(lambda elem: uniprot_common[elem])
    df_final['Common Gene B'] = df_final.loc[:, 'UniProt B'].apply(lambda elem: uniprot_common[elem])
    df_final['Systematic Gene A'] = df_final.loc[:, 'UniProt A'].apply(lambda elem: uniprot_systematic[elem])
    df_final['Systematic Gene B'] = df_final.loc[:, 'UniProt B'].apply(lambda elem: uniprot_systematic[elem])
    print(df_final)
    print('here6')
    df_final = add_surface_residues(pdb_id, df_final)
    print('here7')
    # Writes dataframe to file and returns it for convenience
    write_df(df_final, pdb_id.upper() + "_analysis.csv")
    return df_final