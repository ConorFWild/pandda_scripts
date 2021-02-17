class Constants:
    EXECUTABLE = (
        "#!/bin/bash \n"
        ". /data/share-2/conor/anaconda3/etc/profile.d/conda.sh\n" 
        "conda activate env_rdkit\n"
        "source ~/.bashrc \n"
        "source /data/share-2/conor/xtal_software/ccp4-7.1/bin/ccp4.setup-sh \n"
        "source /data/share-2/conor/xtal_software/phenix/phenix-1.18.2-3874/phenix_env.sh \n"
        "source /data/share-2/conor/xtal_software/buster-2.10/setup.sh \n"
        "python /data/share-2/conor/pandda/pandda_scripts/autobuild/autobuild {model} {xmap} {mtz} {smiles} {x} {y} {z} {out_dir}"
    )

    EXECUTABLE_SCRIPT_FILE = "{dtag}_{event_idx}.sh"

    LOG_FILE = "{event_id}.log"
    OUTPUT_FILE = "{event_id}.out"
    ERROR_FILE = "{event_id}.err"
    REQUEST_MEMORY = "20G"
    JOB = (
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

    JOB_SCRIPT_FILE = "{dtag}_{event_idx}.job"

    COMMAND = "condor_q {job_script_file}"

    MASKED_PDB_FILE = "masked.pdb"

    TRUNCATED_EVENT_MAP_FILE = "truncated.ccp4"

    CUT_EVENT_MAP_FILE = "cut.ccp4"

    LIGAND_PREFIX = "ligand"
    LIGAND_CIF_FILE = "ligand.cif"
    ELBOW_COMMAND = "cd {out_dir}; phenix.elbow {smiles_file} --output=\"{prefix}\"; cd -"

    PANDDA_RHOFIT_SCRIPT_FILE = "/data/share-2/conor/pandda/pandda_scripts/pandda_rhofit.sh"
    RHOFIT_COMMAND = (
        "#!/bin/bash \n"
        "source ~/.bashrc \n"
        ". /data/share-2/conor/anaconda3/etc/profile.d/conda.sh\n" 
        "conda activate env_rdkit\n"  
        "source /data/share-2/conor/xtal_software/ccp4-7.1/bin/ccp4.setup-sh \n"
        "source /data/share-2/conor/xtal_software/phenix/phenix-1.18.2-3874/phenix_env.sh \n"
        "source /data/share-2/conor/xtal_software/buster-2.10/setup.sh \n"
        "{pandda_rhofit} -map {event_map} -mtz {mtz} -pdb {pdb} -cif {cif} -out {out_dir}"
    )
