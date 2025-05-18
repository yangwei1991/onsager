# -*- coding: utf-8 -*-
import os
import itertools
import numpy as np
import mdtraj as mt

def get_atom_indices(traj, resname, element):
    indices = []
    for atom in traj.top.atoms:
        if atom.residue.name == resname and atom.name.startswith(element):
            indices.append(atom.index)
    return indices

def save_structure(name, traj, frame_idx, li_index, all_indices, cutoff=0.2625):
    """
    Save the Li complex structure including full residues of coordinated ligands.

    Parameters:
    - name (str): File name for the saved structure.
    - traj (mdtraj.Trajectory): Trajectory object.
    - frame_idx (int): Frame index to extract the structure.
    - li_index (int): Index of the specific Li atom being analyzed.
    - all_indices (list): All atom indices to consider for distance calculations.
    - cutoff (float): Distance cutoff for including atoms in the complex.
    """
    t = traj[frame_idx]  # Extract the current frame
    pairs = list(itertools.product([li_index], all_indices))
    dists = mt.compute_distances(t, pairs).flatten()
    
    # Find all indices within the cutoff
    indices_within_cutoff = np.argwhere(dists < cutoff).flatten()
    atoms_within_cutoff = [all_indices[i] for i in indices_within_cutoff]

    # Include all atoms from the residues of coordinated atoms
    residues = {atom.residue for atom in traj.top.atoms if atom.index in atoms_within_cutoff}
    full_residue_indices = {li_index}  # Start with the Li atom
    for residue in residues:
        full_residue_indices.update([atom.index for atom in residue.atoms])

    # Slice the trajectory to include only these atoms
    indices = sorted(full_residue_indices)
    t = t.atom_slice(indices)

    # Shift coordinates to the origin
    t.xyz[0, :, :] -= np.min(t.xyz[0, :, :], axis=0)

    # Clean up atom names by removing digits
    for atom in t.top.atoms:
        atom.name = ''.join(filter(lambda x: not x.isdigit(), atom.name))
    
    os.makedirs(os.path.dirname(name), exist_ok=True)
    t.save_xyz(name + '.xyz')

if __name__ == "__main__":
    # Load coordination data
    coordination_data = np.load("coordination_data.npz")
    n1 = coordination_data['n1']
    n2 = coordination_data['n2']
    n3 = coordination_data['n3']
    n4 = coordination_data['n4']

    # Reload trajectory and topology files
    top = '../../entire_box_Additive_MS.prmtop'
    netcdf = '../../entire_box_Additive_MS_nvt_prod2.netcdf'
    traj = mt.load(netcdf, top=top, stride=1)

    # Get atom indices for the respective residues and elements
    LI_indices = get_atom_indices(traj, resname='LI', element='LI')
    DOL_indices = get_atom_indices(traj, resname='DOL', element='O')
    DME_indices = get_atom_indices(traj, resname='DME', element='O')
    TFS_indices = get_atom_indices(traj, resname='TFS', element='O')
    MS1_indices = get_atom_indices(traj, resname='MS1', element='O')

    all_indices = sorted(DOL_indices + DME_indices + TFS_indices + MS1_indices)
    cutoff = 0.2625
    c1, c2, c3, c4 = 0, 0, 0, 0

    # Save structures based on the coordination numbers
    for num, (i, j, k, l) in enumerate(zip(n1, n2, n3, n4)):
        frame_idx = num % traj.n_frames  # Frame index within the trajectory
        li_idx = LI_indices[num // traj.n_frames]  # Li atom index

        if (i, j, k, l) == (0, 6, 0, 0) and c1 < 10:
            c1 += 1
            save_structure(f'MS0600/{c1}', traj, frame_idx, li_idx, all_indices, cutoff)
        if (i, j, k, l) == (0, 4, 0, 1) and c2 < 10:
            c2 += 1
            save_structure(f'MS0401/{c2}', traj, frame_idx, li_idx, all_indices, cutoff)
        if (i, j, k, l) == (0, 5, 0, 0) and c3 < 10:
            c3 += 1
            save_structure(f'MS0500/{c3}', traj, frame_idx, li_idx, all_indices, cutoff)
        if (i, j, k, l) == (0, 2, 0, 2) and c4 < 10:
            c4 += 1
            save_structure(f'MS0202/{c4}', traj, frame_idx, li_idx, all_indices, cutoff)

        if c1 == 10 and c2 == 10 and c3 == 10 and c4 == 10 :
            break
