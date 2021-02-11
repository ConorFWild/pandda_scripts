#!/bin/bash

model=""
xmap=""
mtz=""
smiles=""
x=""
y=""
z=""
out_dir=""

python /data/share-2/conor/pandda/pandda_scripts/batch_cluster.py "$model" "$xmap" "$mtz" "$smiles" "$x" "$y" "$z" "$out_dir"

