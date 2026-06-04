#!/bin/bash

export GMX_MAXBACKUP=-1
# Folder containing PDB files
PDB_FOLDER="."

# Force field to use (change if needed)
FF="amber99sb-ildn"

# Water model (none needed if just fixing PDB)
WATER="none"

# Loop over all PDB files
for pdb in "$PDB_FOLDER"/*.pdb; do
    echo "Processing $pdb ..."

    # Temporary output from pdb2gmx
    temp_output="${pdb%.pdb}_H.pdb"

    # Run pdb2gmx to add hydrogens
    # -f: input PDB, -o: output PDB, -ff: force field, -water: no water
    # Use 'echo 1' to select the first force field if multiple prompts appear
    echo -e "1\n" | gmx pdb2gmx -f "$pdb" -o "$temp_output" -ff "$FF" -water "$WATER" -ignh -missing

    # Overwrite original PDB with hydrogen-added PDB
    mv "$temp_output" "$pdb"

    echo "Hydrogens added and original file overwritten: $pdb"
done

echo "All PDBs processed with hydrogens added."

