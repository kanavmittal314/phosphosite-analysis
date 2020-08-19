# Runs on PyMol to find the surface residues for each PDB ID
# Note: needs to be run like this: pymol -cd "run pymol_surface_resi.py"

from pymol import cmd
from FindSurfaceResidues import findSurfaceResidues
import pandas as pd
import os 
import signal

class TimedOutExc(Exception):
    pass

# Adds a time limit for analyze
def deadline(timeout, *args):
    def decorate(f):
        def handler(signum, frame):
            raise TimedOutExc()

        def new_f(*args):
            signal.signal(signal.SIGALRM, handler)
            signal.alarm(timeout)
            return f(*args)
            signal.alarm(0)

        new_f.__name__ = f.__name__
        return new_f
    return decorate

# Function to analyze the surface residues for a PDB ID
@deadline(3600)
def analyze_surface_residues(pdb_id):
    print('Analyzing surface residues for', pdb_id)
    cmd.fetch(pdb_id)
    findSurfaceResidues(pdb_id, pdb_id=pdb_id, oxygen=True) # Calls on findSurfaceResidues from FindSurfaceResidues.py
    cmd.delete('all')
    
    # Deletes the .cif file created by PyMol
    path = pdb_id.lower() + '.cif'
    if os.path.exists(path):
        os.remove(path)

# Analyzes surface residues for all given PDB IDs and compiles a list of problematic PDB IDs (ones that did not run properly)
def analyze_surface_residues_all(pdb_ids):
    problematic = []
    
    for pdb_id in pdb_ids:
        if not os.path.exists('./surface_residues/' + pdb_id + '_surfaceresidues.csv'):
            try:
                analyze_surface_residues(pdb_id)
            except:
                problematic.append(pdb_id)
    return problematic


if __name__ == '__main__':
    pass
    
    #analyze_surface_residues_all(['6R6G', '4V8Z']) #PDB IDs that did not successfully run within 1 hour
print(analyze_surface_residues_all(pd.read_csv('pdb_ids.csv', header=None, dtype=str).iloc[0].to_list()))


