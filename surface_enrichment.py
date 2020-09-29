# enrichment of phosphosites i

# calculate percent of phosphosites that are on the surface for a given PDB ID
# calculate total percent that are surface_residues
from biopandas.pdb import PandasPdb
import pandas as pd

def get_percent_all_surface(pdb, surface_residues):

    df = pdb[pdb['atom_name'] == 'CA']
    return len(surface_residues)/len(df)

def get_percent_phosphosites_surface(phos_analysis):
    
    return len(set(zip(phos_analysis[phos_analysis['Surface A']]['Chain A'], phos_analysis[phos_analysis['Surface A']]['Residue A'])) | set(zip(phos_analysis[phos_analysis['Surface B']]['Chain B'], phos_analysis[phos_analysis['Surface B']]['Residue B']))) / len(set(zip(phos_analysis['Chain A'], phos_analysis['Residue A'])) | set(zip(phos_analysis['Chain B'], phos_analysis['Residue B'])))

def calculate_enrichment(pdb_id):
    all = get_percent_all_surface(PandasPdb().fetch_pdb(pdb_id).df['ATOM'], pd.read_csv('./surface_residues/' + pdb_id + '_surfaceresidues.csv'))
    phos = get_percent_phosphosites_surface(pd.read_csv('./phosphosite_analysis_results/' + pdb_id + '_analysis.csv'))
    print(all)
    print(phos)
    print(phos/all)
    return phos/all

def calculate_enrichments(pdb_ids):
    for pdb_id in pdb_ids:
        print(calculate_enrichment(pdb_id))


print(calculate_enrichments(['1M2O', '1M38', '1M0T']))

#must anyways re-run phosphosite_analysis.py as the surface_residues are the pre-oxygen ones