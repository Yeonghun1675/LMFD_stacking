import networkx as nx
import numpy as np
from itertools import chain
import ase
import ase.io
import ase.data
import ase.visualize
import ase.neighborlist
try:
    from ase.utils import natural_cutoffs
except:
    from ase.neighborlist import natural_cutoffs

from scipy.sparse.csgraph import connected_components


def covalent_neighbor_list(atoms, scale=1.2, neglected_species=(), neglected_indices=()):
    cutoffs = natural_cutoffs(atoms, mult=scale)
    
    return ase.neighborlist.neighbor_list('ijD', atoms, cutoff=cutoffs)


def get_connected_component(atoms):
    I, J, _ = covalent_neighbor_list(atoms)
    graph = nx.Graph(zip(I, J))

    ccs = list((nx.connected_components(graph)))
    return ccs


def get_max_height(mol, cutoff=20, height=60):
    atoms = ase.io.read(mol)
    con_comps = get_connected_component(atoms)
    component = chain.from_iterable([cc for cc in con_comps if len(cc) > cutoff])
    new_atoms = atoms[list(component)]
    z_position = new_atoms.get_positions()[:,2]
    z_position = z_position[z_position< height]
    max_height = np.max(z_position)
    
    return max_height


if __name__ == '__main__':
    mol = sys.argv[1]
    max_height = get_max_height(mol)
    print (max_height)
