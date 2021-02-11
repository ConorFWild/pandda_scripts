class Constants:
    EXECUTABLE = (
        ""
    )

    EXECUTABLE_SCRIPT_FILE = "{dtag}_{event_idx}.sh"

    JOB = (
        ""
    )

    JOB_SCRIPT_FILE = "{dtag}_{event_idx}.job"

    COMMAND = "condor_q {job_script_file}"

    MASKED_PDB_FILE = "masked.pdb"

    TRUNCATED_EVENT_MAP_FILE = "truncated.ccp4"

    LIGAND_PREFIX = "ligand"
    LIGAND_CIF_FILE = "ligand.cif"
    ELBOW_COMMAND = "cd {out_dir}; phenix. {smiles_file} --output=\"{autobuilding_ligand}\""

    PANDDA_RHOFIT_SCRIPT_FILE = "/data/share-2/conor/pandda/pandda_scripts/pandda_rhofit.sh"
    RHOFIT_COMMAND = (
        "#!/bin/bash \n"
        "conda activate pandda \n"
        "source ~/.bashrc \n"
        "source /data/share-2/conor/xtal_software/ccp4-7.1/bin/ccp4.setup-sh \n"
        "source /data/share-2/conor/xtal_software/phenix/phenix-1.18.2-3874/phenix_env.sh \n"
        "source /data/share-2/conor/xtal_software/buster-2.10/setup.sh \n"
        "{pandda_rhofit} -map {event_map} -mtz {mtz} -pdb {pdb} -cif {cif} -out {out_dir_path}"
    )
