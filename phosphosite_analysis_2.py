from biopandas.pdb import PandasPdb
import math
import pandas as pd
import numpy as np
from scipy.spatial.distance import pdist
import itertools
import argparse
import os
from io import StringIO

#Set up command-line arguments
def get_args():
  parser = argparse.ArgumentParser(description='Run phosphosite distance analysis')
  parser.add_argument('phosphosites', type = str)
  parser.add_argument('--path', type = str)
  parser.add_argument('--rcsb', type = str)
  return parser.parse_args()

#Uses four-letter pdb code to get .pdb file from rcsb.org. Returns an error if invalid code.
def get_rcsb(pdb):
  print("RCSB...")
  try:
    ppdb = PandasPdb().fetch_pdb(pdb)
    print("Fetched")
    return ppdb
  except:
    print("PDB file not found on RCSB")
    raise SystemExit

#Gets the PDB file from local user 
def get_file_input(path):
  try:
    return PandasPdb().read_pdb(path)
  except:
    print("No PDB file found at " + path)
    raise SystemExit

#Determines where to get pdb file from RCSB or local drive
def get_ppdb(args):
  if args.path:
    return get_file_input(args.path)
  elif args.rcsb:
    return get_rcsb(args.rcsb)
  else:
    print("No input file or PDB code")
    raise SystemExit

#Gets phosphosite list from local drive
def get_phosphosites(args):
  path = args.phosphosites
  try:
    df = pd.read_csv(path, squeeze = True)
    print(df)
    return df
  except:
    print("No CSV file found for phosphosites")
    raise SystemExit

#Filters the atoms to keep only the atoms that are on the appropriate residues and are the correct oxygens
def get_phospho_coords(residue_nums, df):
  residue_criteria = df['residue_number'].isin(residue_nums)
  serine_criteria = (df['residue_name'] == 'SER') & (df['atom_name'] == 'OG')
  threonine_criteria = (df['residue_name'] == 'THR') & (df['atom_name'] == 'OG1')
  tyrosine_criteria = (df['residue_name'] == 'TYR') & (df['atom_name'] == 'OH')
  alt_loc_criteria = (df['alt_loc'] == '') | (df['alt_loc'] == 'A')
  #return df[serine_criteria | threonine_criteria | tyrosine_criteria]
  return df[residue_criteria & (serine_criteria | threonine_criteria | tyrosine_criteria) & alt_loc_criteria]

'''
first get the residue_pairs from all possible phosphosites and calculate linear distance
then at each residue_pairs, calculate euclidean distance; 
  still create the phospho_coords df and search through there, instead of the entire df
chains should be easily accounted for in this method, as will be missing residues

still to be figured out: discrepancies (use SEQRES and SEQADV in PDB file)
'''
#Gives pairwise combinations of all residues
def get_residue_pairs(df):
  return np.asarray(list(itertools.combinations(df['residue_number'].tolist(), 2))), np.asarray(list(itertools.combinations(df['chain_id'].tolist(), 2)))
  # return "", np.asarray(list...df['chain_id'].tolist(), 2)

#Gives pairwise euclidean distances of all residues
def get_euclidean_distances(df):
  coords = df.iloc[:,11:14]
  return pdist(coords.to_numpy())

def get_linear_distance():
  pass



#Writes the final dataframe to an output file
def write_df(df, file_name):
    dir = "./phosphosite_analysis_results"
    path = dir + "/" + file_name
    if not os.path.exists(dir):
      os.mkdir(dir)
    if os.path.exists(path):
      os.remove(path)
    df.to_csv(path, index = False)

def get_input():
  args = get_args()
  ppdb = get_ppdb(args)
  phosphosites = get_phosphosites(args)
  return ppdb, phosphosites

'''
def phosphosites_missing(df, phosphosites):
  return [elem for elem in phosphosites if elem not in df['residue_number'].values]
'''

def gene_chain(df)

def split_by_chain(df):
  return df.groupby(['chain_id'])

def glinear_distance(dff, phosphosites, chain):

  #print(itertools.combinations(df['residue_number'].tolist()))
  #print(np.asarray(list(itertools.combinations(df['residue_number'].tolist(), 2))))
  y = np.hsplit(np.asarray(list(itertools.combinations(phosphosites, 2))), 2)
  #print("Y")
  print(y)
  dff = pd.DataFrame.from_dict({"Residue A": np.squeeze(y[0]), "Residue B": np.squeeze(y[1])})#, "Linear Distance": np.squeeze(np.absolute(np.subtract(y[0], y[1])))})
  dff.insert(0, 'Chain A', chain)
  dff.insert(2, 'Chain B', chain)
  return dff

def process_pdb(ppdb):


#Main part of program
def analyze_phosphosite_distances(ppdb, phosphosites):
  #Gets and filters the pdb file to include only relevant atoms
  process_ppdb(ppdb)

  phosphosites.sort()
  df = ppdb.df['ATOM']
  
  df_filtered = get_phospho_coords(phosphosites, df)
  missing_phosphosites = phosphosites_missing(df_filtered, phosphosites)
  print(missing_phosphosites)
  included_phosphosites = df_filtered['residue_number'].values
  print(included_phosphosites)

  xx = (split_by_chain(df))
  print("here")
  
  dff = pd.DataFrame(columns = ['Chain A', 'Residue A', 'Chain B', 'Residue B'])
  print(xx.groups.keys())
  for i in xx.groups.keys():
    print(xx.get_group(i))
    dff = pd.concat([dff, glinear_distance(xx.get_group(i), phosphosites, i)])
  print(dff)

  #list 1: missing, list 2: all, all combintions for all possible chains?, add in linear distances!
  #methionine re-alignment

  #Gets all pairwise combinations and splits them into two columns for residue_1 and residue_2
  residue_pairs, chain_pairs = get_residue_pairs(df_filtered)
  residue_pairs_split, chain_pairs_split = np.hsplit(residue_pairs, 2), np.hsplit(chain_pairs, 2)
  residue_1, chain_1 = residue_pairs_split[0], chain_pairs_split[0]
  residue_2, chain_2 = residue_pairs_split[1], chain_pairs_split[1]

  #Gets euclidean and linear distances
  euclidean_distances = get_euclidean_distances(df_filtered)
  linear_distances = np.absolute(np.subtract(residue_1, residue_2))

  #Creates a final dataframe with the two residues, 3D euclidean distance, and linear distance
  dictionary = {"Chain A": np.squeeze(chain_1), "Residue A": np.squeeze(residue_1), "Chain B": np.squeeze(chain_2), "Residue B": np.squeeze(residue_2), "3D distance": euclidean_distances}
  df_final = pd.DataFrame.from_dict(dictionary)

  #df_final.insert(4, "Linear dISTANCES", None)

  print("Merging")
  df_final = df_final.merge(dff, how='outer', left_on=['Chain A', 'Residue A', 'Chain B', 'Residue B'], right_on = ['Chain A', 'Residue A', 'Chain B', 'Residue B'])
  #print(df_final.merge(dff, how='outer', left_on=['Chain A', 'Residue A', 'Chain B', 'Residue B'], right_on = ['Chain B', 'Residue B', 'Chain A', 'Residue A']))

  df_final['Linear Distances'] = np.absolute(np.subtract(df_final['Residue A'], df_final['Residue B']))

  #Writes dataframe to file and returns it for convenience
  write_df(df_final, "_analysis.csv")
  return df_final

ppdb, phosphosites = get_input()
print(analyze_phosphosite_distances(ppdb, [43, 719, 24, 21, 223, 225]))




