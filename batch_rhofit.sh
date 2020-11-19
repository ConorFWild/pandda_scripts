#!/bin/bash 

# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/data/share-2/conor/anaconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/data/share-2/conor/anaconda3/etc/profile.d/conda.sh" ]; then
        . "/data/share-2/conor/anaconda3/etc/profile.d/conda.sh"
    else
        export PATH="/data/share-2/conor/anaconda3/bin:$PATH"
    fi
fi
unset __conda_setup


conda activate pandda 
source ~/.bashrc 
source /data/share-2/conor/xtal_software/ccp4-7.1/bin/ccp4.setup-sh 
source /data/share-2/conor/xtal_software/phenix/phenix-1.18.2-3874/phenix_env.sh 
source /data/share-2/conor/xtal_software/buster-2.10/setup.sh 
python /data/share-2/conor/pandda/pandda_scripts/batch_rhofit.py --pandda_dirs_dir=/data/share-2/conor/pandda/output/test_batch_pandda/ --autobuild_dirs_dir=/data/share-2/conor/pandda/output/autobuild