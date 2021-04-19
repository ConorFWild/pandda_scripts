#!/bin/bash
cutoff=$1
root="/data/share-2/conor/pandda"
data_dirs="$root/data/pandda_inputs/XX02KALRNA_$cutoff"
output_dir="$root/output/test_local_cluster_XX02KALRNA_$cutoff"
submit_script="$output_dir/test_local_cluster_XX02KALRNA_$cutoff.sh"
job_script="$output_dir/test_local_cluster_XX02KALRNA_$cutoff.job"
output_file="$output_dir/test_local_cluster_XX02KALRNA_$cutoff.out"
error_file="$output_dir/test_local_cluster_XX02KALRNA_$cutoff.err"
log_file="$output_dir/test_local_cluster_XX02KALRNA_$cutoff.log"
local_pandda_cluster_script="$root/pandda_local_cluster/pandda_local_cluster.py"


mkdir $output_dir

echo -e "#!/bin/bash \n/data/share-2/conor/anaconda3/envs/lab/bin/python $local_pandda_cluster_script $data_dirs $output_dir --cutoff=$cutoff" > $submit_script
echo -e "# Executable script at $submit_script"

chmod 777 $submit_script

echo -e "Executable = $submit_script \noutput = $output_file \nerror = $error_file \nlog = $log_file \nrequest_memory = 100 GB \nQueue" > $job_script
echo -e "# Job script at $job_script"

condor_submit $job_script
