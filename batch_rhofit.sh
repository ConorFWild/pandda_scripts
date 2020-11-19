#!/bin/bash \n
conda activate pandda 
source ~/.bashrc 
source /data/share-2/conor/xtal_software/ccp4-7.1/bin/ccp4.setup-sh 
source /data/share-2/conor/xtal_software/phenix/phenix-1.18.2-3874/phenix_env.sh 
source /data/share-2/conor/xtal_software/buster-2.10/setup.sh 
python /data/share-2/conor/pandda/pandda_scripts/batch_rhofit.py