class Constants:
    SKELETON_SCORE_FILE = "score.json"

    SCORE_SCRIPT = (
        "#!/bin/bash \n"
        ". /data/share-2/conor/anaconda3/etc/profile.d/conda.sh\n"
        "conda activate lab\n"
        "source ~/.bashrc \n"
        "source /data/share-2/conor/xtal_software/ccp4-7.1/bin/ccp4.setup-sh \n"
        "source /data/share-2/conor/xtal_software/phenix/phenix-1.18.2-3874/phenix_env.sh \n"
        "source /data/share-2/conor/xtal_software/buster-2.10/setup.sh \n"
        "/data/share-2/conor/anaconda3/envs/lab/bin/python /data/share-2/conor/pandda/pandda_scripts/skeleton/skeleton.py {structure} {event_map} {out_dir}"
    )

    LOG_FILE = "{build_id}.log"
    OUTPUT_FILE = "{build_id}.out"
    ERROR_FILE = "{build_id}.err"
    JOB_SCRIPT = (
        "#################### \n"
        "# \n"
        "# Example 1                                   \n"
        "# Simple HTCondor submit description file \n"
        "#                          \n"
        "####################    \n"

        "Executable   = {executable_file} \n"
        "Log          = {log_file} \n"
        "Output = {output_file} \n"
        "Error = {error_file} \n"

        "request_memory = {request_memory} GB \n"

        "Queue"
    )

    SCORE_SCRIPT_FILE = "{build_id}.sh"
    JOB_SCRIPT_FILE = "{build_id}.job"

    SUBMIT_COMMAND = "condor_submit {job_script_file}"

    REQUEST_MEMORY = "20"

    DEBUG = True

