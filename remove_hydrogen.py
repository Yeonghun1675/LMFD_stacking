import sys
import random
import ase
import ase.io
from pathlib import Path

def remove_hydrogen(molecule, output='active.xyz', num=4):
    atoms = ase.io.read(f'{molecule}.xyz')
    
    if num:
        hydrogens = [atom.index for atom in atoms if atom.symbol == 'H']
        delete_list = random.sample(hydrogens, num)
        del atoms[delete_list]
    
    ase.io.write(output, atoms)

if __name__ == '__main__':
    molecule = sys.argv[1]
    remove_hydrogen(molecule)


