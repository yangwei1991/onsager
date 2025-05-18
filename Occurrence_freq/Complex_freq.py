# -*- coding: utf-8 -*-

import itertools
import numpy as np
import mdtraj as mt

def get_atom_indices(traj, resname, element):
    indices = []
    for atom in traj.top.atoms:
        if atom.residue.name == resname and atom.name.startswith(element):
            indices.append(atom.index)
    return indices

# Load the trajectory and topology files
top = '../../entire_box_Additive_MS.prmtop'
netcdf = '../../entire_box_Additive_MS_nvt_prod2.netcdf'
traj = mt.load(netcdf, top=top, stride=1)

# Get atom indices for the respective residues and elements
LI_indices = get_atom_indices(traj, resname='LI', element='LI')
DOL_indices = get_atom_indices(traj, resname='DOL', element='O')
DME_indices = get_atom_indices(traj, resname='DME', element='O')
TFS_indices = get_atom_indices(traj, resname='TFS', element='O')
MS1_indices = get_atom_indices(traj, resname='MS1', element='O')

# Define the cutoff distance (in nm)
cutoff = 0.2625

# Initialize lists to store counts for each ligand
n1, n2, n3, n4 = [], [], [], []

# Loop over each Li ion and count ligands within the cutoff distance
for ndx in LI_indices:
    # DOL ligands around Li
    dists = mt.compute_distances(traj, itertools.product([ndx], DOL_indices))
    n1 += list(np.sum(dists < cutoff, axis=1))

    # DME ligands around Li
    dists = mt.compute_distances(traj, itertools.product([ndx], DME_indices))
    n2 += list(np.sum(dists < cutoff, axis=1))

    # TFS ligands around Li
    dists = mt.compute_distances(traj, itertools.product([ndx], TFS_indices))
    n3 += list(np.sum(dists < cutoff, axis=1))

    # MS1 ligands around Li
    dists = mt.compute_distances(traj, itertools.product([ndx], MS1_indices))
    n4 += list(np.sum(dists < cutoff, axis=1))

# Create a dictionary to store occurrences of each complex type [n1, n2, n3, n4]
count_dict = {}
for i, j, k, l in zip(n1, n2, n3, n4):
    key = (i, j, k, l)  # Create a tuple as the key
    count_dict[key] = count_dict.get(key, 0) + 1  # Increment count

# Normalize the counts based on the number of frames and Li ions
vals = np.array(list(count_dict.values())) / (traj.n_frames * len(LI_indices))

# Sort the complexes by their occurrence frequency
indices = np.flip(np.argsort(vals))
keys = np.array(list(count_dict.keys()))

# Print the occurrence of each complex in the format [n1, n2, n3, n4] occurrence
for k, v in zip(keys[indices], vals[indices]):
    print(f"{k}: {v:.4f}")

# Save coordination numbers to a file
np.savez("coordination_data.npz", n1=n1, n2=n2, n3=n3, n4=n4)