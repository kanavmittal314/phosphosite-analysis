# Input = list of all PDB IDs
# Fetch all PDBs and feed into analyze_phosphosite_distances
# Get all chain_genes for each ppdb; dictionary mapping each PDB ID to a chain_gene
# Go through dict, at each entry of chain_gene, get a list of phosphosites for that UniProt ID
# Apply operation to phosphosites to make it so that is in from chain: phosphosites using chain_gene_dict
# Create a phosphosite_dict from those phosphosites and feed into analyze_phosphosite_distances

from phosphosite_analysis_7_16 import get_chain_gene, analyze_phosphosite_distances, get_rcsb, get_phosphosites, get_chain_gene_dict, get_fasta_dict
import pandas as pd
import urllib.parse
import urllib.request

# Running list of all PDB IDs considered to be problematic
problematic = []

# Gets the ppdbs and the chain_gene_dicts for all pdb_ids
def get_info(pdb_ids):
    ppdb_dict = {}
    chain_gene_dict = {}

    for i in pdb_ids:
        ppdb = get_rcsb(i)
        ppdb_dict[i] = ppdb
        chain_gene_dict[i] = get_chain_gene_dict(ppdb)
    
    return ppdb_dict, chain_gene_dict

# Gets the phosphosites for each PDB ID
def get_all_phos(pdb_ids, chain_gene_dict):

    total_phos = {}
    global problematic
    for i in pdb_ids:
        chain_gene = chain_gene_dict[i]
        phosphosites = {}
        try:
            for j in chain_gene:
                phosphosites[j] = get_phosphosites("phosphosites/" + chain_gene[j] + "_phos.csv").to_numpy()
        except SystemExit:
            print("Problem with getting phosphosites for PDB ID", i)
            problematic.append(i)
        total_phos[i] = phosphosites
    
    return total_phos

# Runs the phosphosite_analysis on each PDB using the phosphosites from get_all_phos()
def analyze_all(pdb_ids, ppdb_dict, total_phos, uniprot_common_mapping, uniprot_systematic_mapping, sgd_mapping, fasta_dict):
    global problematic
    for i in pdb_ids:
        try:
            print(analyze_phosphosite_distances(i, ppdb_dict[i], total_phos[i], uniprot_common_mapping, uniprot_systematic_mapping, sgd_mapping, fasta_dict))
        except:
            print("Problem with analyzing PDB ID", i)
            problematic.append(i)

def get_sgd_mapping(uniprot_ids):
    url = 'https://www.uniprot.org/uploadlists/' 
    
    params = {
    'from': 'ACC+ID',
    'to': 'SGD_ID',
    'format': 'tab',
    'query': ' '.join(uniprot_ids)
    }

    data = urllib.parse.urlencode(params)
    data = data.encode('utf-8')
    req = urllib.request.Request(url, data)
    with urllib.request.urlopen(req) as f:
        response = f.read()
    txt = (response.decode('utf-8').split())[2:]
    print(len(txt))
    a = {}
    for i in range(len(txt)):
        if not txt[i].startswith('S') and txt[i+1].startswith('S'):
            a[txt[i]] = txt[i+1]
    print(a)
    print(len(a.keys()))
    
    
    #raise SystemExit

    
    return a
    

# Main function that calls all other functions and writes the problematic PDB IDs to a file
def solve(pdb_ids):

    #sgd_ids, fasta_dict, Lanz_common_systematic.csv ==> args of analyze_phosphosite_distances
    #calculated in yeast_proteins.solve()

    common_systematic_gene_names = pd.read_csv("Lanz_systematic_common_gene.csv")
    common_systematic_gene_names.loc[:, 'Uniprot_id'].to_csv("All_UniProt_ID.csv", header=False, index=False)
    
    #print(common_systematic_gene_names.to_dict())

    uniprot_common_mapping = dict(zip(common_systematic_gene_names.loc[:,'Uniprot_id'], common_systematic_gene_names.loc[:,'Gene(s)']))
    #print(uniprot_common_mapping)

    uniprot_systematic_mapping = dict(zip(common_systematic_gene_names.loc[:,'Uniprot_id'], common_systematic_gene_names.loc[:,'SystematicGeneName']))
    #print(uniprot_systematic_mapping)
    
    #print(common_systematic_gene_names.loc[:,'Uniprot_id'])
    #print(len(common_systematic_gene_names.loc[:,'Uniprot_id']))
    #print(len(set(common_systematic_gene_names.loc[:,'Uniprot_id'])))

    

    sgd_mapping = get_sgd_mapping(common_systematic_gene_names.loc[:,'Uniprot_id'])

    
    print(sgd_mapping['Q08220'])
    print(sgd_mapping)
    
    ppdb_dict, chain_gene_dict = get_info(pdb_ids)

    total_phos = get_all_phos(pdb_ids, chain_gene_dict)

    fasta_dict = get_fasta_dict()

    analyze_all(pdb_ids, ppdb_dict, total_phos, uniprot_common_mapping, uniprot_systematic_mapping, sgd_mapping, fasta_dict)

    pd.Series(problematic).to_csv("problematic_pdb_ids.csv")
    return "Finished"




#print(solve(['1m0t']))
#raise SystemExit
print(solve(['1m0t', '1m2o', '1m38']))
#print(solve(['1m2o', '1m38']))
