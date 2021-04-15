#!/bin/bash
root="/data/share-2/conor/pandda"
data_dirs="$root/data/pandda_inputs/HAO1A"
output_dir="$root/output/test_local_pandda_HAO1A"
submit_script="$output_dir/test_local_pandda.sh"
job_script="$output_dir/test_local_pandda.job"
output_file="$output_dir/HAO1A.out"
error_file="$output_dir/HAO1A.err"
log_file="$output_dir/HAO1A.log"
local_pandda_script="$root/local_pandda/local_pandda_ligand_gpu.py"
phenix_setup="/data/share-2/conor/xtal_software/phenix/phenix-1.18.2-3874/phenix_env.sh"
rhofit_setup="/data/share-2/conor/xtal_software/buster-2.10/setup.sh"


mkdir $output_dir

echo -e "#!/bin/bash \ncd $root; /data/share-2/conor/anaconda3/envs/gpu/bin/python $local_pandda_script $data_dirs $output_dir" > $submit_script

chmod 777 $submit_script

echo -e "Executable = $submit_script \noutput = $output_file \nerror = $error_file \nlog = $log_file \nrequest_memory = 10 GB \nQueue" > $job_script

condor_submit $job_script
