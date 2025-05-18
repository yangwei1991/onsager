#!/bin/bash

# Loop through Conformer directories
for a in {3..20}; do
    cd "Conformer${a}" || exit
    
    # Loop through PDB files
    for i in $(seq -w 1 27); do
        # Create tleap input file
        cat << EOF > "tleap${i}.in"
source leaprc.protein.ff14SB
source leaprc.gaff2
source leaprc.water.tip3p
TFS = loadmol2 TFSI.mol2
loadamberparams TFSI.frcmod
loadamberparams frcmod.ions234lm_126_tip3p
mol = loadpdb s0000${i}.pdb
saveamberparm mol s0000${i}_dry.prmtop s0000${i}_dry.inpcrd
quit
EOF
    tleap -s -f tleap${i}.in > tleap${i}.out
    done
    cd ..
done




