#!/bin/bash

source ~/.bashrc
conda activate pandda

pandda="/data/share-2/conor/pandda/code/pandda_gemmi/pandda_gemmi/analyse.py"
data_dirs="/data/share-2/conor/pandda/data/pandda_inputs/PHIPA"
cluster_cutoff_distance_multiplier=1.1 
min_blob_volume=9.0
out_dir="/data/share-2/conor/pandda/output/test/test_pandda"
pdb_regex="final.pdb"
mtz_regex="final.mtz"


python $pandda --data_dirs=$data_dirs --cluster_cutoff_distance_multiplier=$cluster_cutoff_distance_multiplier --min_blob_volume=$min_blob_volume --out_dir=$out_dir --pdb_regex=$pdb_regex --mtz_regex=$mtz_regex



