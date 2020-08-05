# Runs on PyMol to find the surface residues for each PDB ID
# Note: needs to be run like this: pymol -cd "run pymol_surface_resi.py"

from pymol import cmd
from FindSurfaceResidues import findSurfaceResidues
import pandas as pd
import os 
import signal

class TimedOutExc(Exception):
    pass

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

@deadline(3600)
def analyze(i):
    print(i)
    cmd.fetch(i)
    print('here')
    findSurfaceResidues(i, pdb_id=i)
    print('here2')
    cmd.delete('all')
    print('here3')
    
    path = i.lower() + '.cif'
    if os.path.exists(path):
        os.remove(path)

def solve(pdb_ids):
    problematic = []
    for i in pdb_ids:
        if not os.path.exists('surface_residues/' + i + '_surfaceresidues.csv'):
            try:
                analyze(i)
            except:
                problematic.append(i)
            
    print(problematic)



#solve(pd.read_csv('pdb_ids.csv', header=None, dtype=str).iloc[0].to_list())
solve(['6R6G', '4V8Z']) #PDB IDs that did not successfully run within 1 hour


