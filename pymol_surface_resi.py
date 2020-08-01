# Runs on PyMol to find the surface residues for each PDB ID
# Note: needs to be run like this: pymol -d "run pymol_surface_resi.py"

from pymol import cmd
from FindSurfaceResidues import findSurfaceResidues
import pandas as pd
import os 


def solve(pdb_ids):
    problematic = []
    for i in pdb_ids:
        if not os.path.exists('surface_residues/' + i + '_surfaceresidues.csv'):
            try:
                print(i)
                cmd.fetch(i)
                findSurfaceResidues(i, pdb_id=i)
                cmd.delete('all')
                path = i.lower() + '.cif'
                if os.path.exists(path):
                    os.remove(path)
            except:
                problematic.append(i)
    print(problematic)

solve(pd.read_csv('pdb_ids.csv', header=None, dtype=str).iloc[0].to_list())



