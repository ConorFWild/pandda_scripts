#!/bin/bash 

conda activate pandda 
source ~/.bashrc 
source /data/share-2/conor/xtal_software/ccp4-7.1/bin/ccp4.setup-sh 
source /data/share-2/conor/xtal_software/phenix/phenix-1.18.2-3874/phenix_env.sh 
source /data/share-2/conor/xtal_software/buster-2.10/setup.sh 
/data/share-2/conor/anaconda3/envs/pandda/bin/python /data/share-2/conor/pandda/pandda_scripts/batch_rhofit.py --pandda_dirs_dir=/data/share-2/conor/pandda/output/test_batch_pandda/ --autobuild_dirs_dir=/data/share-2/conor/pandda/output/autobuild