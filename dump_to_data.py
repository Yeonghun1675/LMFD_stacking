#!/usr/bin/env python
# coding: utf-8

import os
import sys
import argparse
import re
import subprocess
import random
import shutil
from itertools import zip_longest
import numpy as np
from remove_hydrogen import remove_hydrogen
from get_connected_atoms import get_max_height


def data_to_info(index):
    data_path = f'./write_data/data.{index}'
    data = ''
    atoms = []
    velocity = []

    with open(data_path) as f:
        for line in f:
            if line.startswith('Atoms # charge'):
                data += line + "\n"
                break
            else:
                data += line
        next(f)
        
        for line in f:
            if line.startswith('Velocities'):
                break
            else:
                blocks = line.strip().split()
                if len(blocks) == 9:
                    atoms.append(blocks)
        next(f)
        
        for line in f:
            blocks = line.strip().split()
            if len(blocks) == 4:
                velocity.append(blocks)

    atoms.sort(key=lambda t: int(t[0]))
    velocity.sort(key=lambda t: int(t[0]))  
    assert len(atoms) == len(velocity)
    assert len(atoms) == int(atoms[-1][0]), f'len atoms : {len(atoms)}, last atoms : {int(atoms[-1][0])}'

    return data, atoms, velocity


def change_atoms_from_data(data, num):
    return re.sub(r"[0-9]+(?= atoms)", str(num), data)


def remove_hydrogen_from_data(atoms, velocity, num_remove=4, sort_by_height=True):
    if not num_remove:
        return atoms, velocity

    if not sort_by_height:
        print ('sort_by_random\n')
        index_hydrogen = [atom[0] for atom in atoms if atom[1] == '2']
        try:
            rm_idx = random.sample(index_hydrogen, num_remove)
        except ValueError:
            return atoms, velocity

        #print (rm_idx)

    else:
        print ('sort_by_height\n')
        index_hydrogen = [(atom[0], atom[5]) for atom in atoms if atom[1] == '2']
        s_ls = sorted(index_hydrogen, key=lambda t: float(t[1]))

        rm_idx = [s[0] for s in s_ls[:num_remove]]
    
    atoms = [atom for atom in atoms if atom[0] not in rm_idx]
    velocity = [v for v in velocity if v[0] not in rm_idx]

    for i, (atom, v) in enumerate(zip(atoms, velocity), 1):
        assert atom[0] == v[0]
        atoms[i-1][0] = str(i)
        velocity[i-1][0] = str(i)
        
    return atoms, velocity
             

def write_xyz_from_data(atoms, output_path):
    atom_index = {'1':'C', '2':'H', '3':'O', '4':'N'}

    with open(output_path, 'w') as f:
        f.write(f"{len(atoms)}\n\n")
        for atom in atoms:
            pos = atom[3:6]
            type_ = atom_index.get(atom[1])
            xyz = [type_] + pos

            f.write("        ".join(xyz) + "\n")


def get_box_from_data(data):
    box = []
    for line in data.split('\n'):
        if re.search(r'[xyz]lo [xyz]hi', line):
            box.append(line.split()[:2])
    assert len(box) == 3
    return box


def xyz_to_data(xyz_file, atoms, velocity, data, output = 'tmp.data'):
    atom_index = {'C':'1', 'H':'2', 'O':'3', 'N':'4'}

    with open(xyz_file) as f:
        num = int(next(f).strip())       # first_line : num atoms
        next(f)                          # pass 1 line
        
        new_atoms = []
        for i, line in enumerate(f):
            line = line.split()
            type_ = atom_index.get(line[0])
            xyz = line[1:]
            vec = [str(i+1), type_, '0.0'] + xyz + ['0.0', '0.0', '0.0']
            assert len(vec) == 9
            new_atoms.append(vec)

    assert num == len(new_atoms)
    
    data = change_atoms_from_data(data, num)

    with open(output, 'w') as f:
        f.write(data) # write Num atoms, types, boxs, and Masses

        for I, J in zip_longest(atoms, new_atoms):
            if I:
                f.write(" ".join(I) + "\n")
            else:
                f.write(" ".join(J) + "\n")

        f.write("\nVelocities\n\n")

        count = 0
        for i in range(num):
            try:
                line = " ".join(velocity[i])
            except IndexError:
                line = f"{i+1} 0.0                    0.0                    0.0"
                if not count:
                    count = i+1

            f.write(line + "\n")
    return count

def make_inp(box, index, LMFD, output='tmp.inp', zlim=None):
    box = np.array(box, dtype='float')

    pipe = [f"sed 's/(output)/{index}/' sample.inp", 
            f"sed 's/(index)/{index}/'",
            f"sed 's/(LMFD)/{LMFD}/'"]
    
    for cd, value in zip(['x', 'y'], box):
        pipe.append(f"sed 's/({cd}1)/{value[0]}/'")
        pipe.append(f"sed 's/({cd}2)/{value[1]-0.5}/'")
    
    if zlim:
        pipe.append(f"sed 's/(z1)/{zlim}/'")
        pipe.append(f"sed 's/(z2)/{zlim+5}/'")
    else:
        pipe.append(f"sed 's/(z1)/{box[2][0]}/'")
        pipe.append(f"sed 's/(z2)/{box[2][1]-0.5}/'")


    pipeline = " | ".join(pipe)
    pipeline += f' > {output}'
                    
    assert subprocess.check_call(pipeline, shell=True) == 0
    

def make_in(index, num, output, eq_state=False):

    if eq_state:
        velocity = ' '
    else:
        velocity = 'velocity        new set NULL NULL -0.002'

    pipe = [f"sed 's/(index)/{index}/' sample.in", 
            f"sed 's/(index+1)/{index+1}/'",
            f"sed 's/(newindex)/{num}/'",
            f"sed 's/(velocity)/{velocity}/'"
            ]
    pipeline = " | ".join(pipe)
    pipeline += f' > {output}'
    assert subprocess.check_call(pipeline, shell=True) == 0



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='dump_to_datafile')
    parser.add_argument('--index', type=int, 
                        help='index')
    parser.add_argument('--num-active', type=int, default=0,
                        help='number_of_active_site')
    parser.add_argument('--molecule', type=str, default='LMFD')
    parser.add_argument('--remove-hydrogen', type=int, default=0)
    parser.add_argument('--no-packmol', action='store_true',)
    parser.add_argument('--sort-by-height', action='store_true',)
    parser.add_argument('--eq-state', action='store_true',)


    args = parser.parse_args()    
    index = args.index
    LMFD = args.molecule
    num_active = args.num_active
    rm_h = args.remove_hydrogen
    
    # GET DATA
    data, atoms, velocity = data_to_info(index)
    box = get_box_from_data(data)

    print (args.sort_by_height, type(args.sort_by_height))
    if rm_h:
        atoms, velocity = remove_hydrogen_from_data(atoms, velocity, num_remove=rm_h, 
        sort_by_height=args.sort_by_height)
    
    # PACKMOL
    tmp_structure = f'tmp_structure/{index}.xyz'
    write_xyz_from_data(atoms, tmp_structure)

    mol_name = LMFD.replace('/', '_')
    remove_hydrogen(f'../molecule/{LMFD}', f'active_{mol_name}.xyz', num=num_active)

    if not args.no_packmol:
        make_inp(box, index, f'active_{mol_name}', output=f"inp/{index}.inp", zlim = 55) 
        with open(f'./inp/{index}.inp', 'rb') as f:
            subprocess.check_call('packmol', stdin=f, shell=True)
    else:
        print ("packmol does not used in these system")
        shutil.copy(tmp_structure, f'./packmol_structure/{index}.xyz')

    
    # MAKE LAMMPS INPUT
    num = xyz_to_data(f'packmol_structure/{index}.xyz', atoms, velocity, data, output=f'data/{index}.data')
    make_in(index, num, output=f"lammps/{index}.in", eq_state=args.eq_state)
