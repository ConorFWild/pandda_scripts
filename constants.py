
class Constants:
    DEBUG: int = 0
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

    JOB_SCRIPT_FILE = "{system_name}.job"

    LOG_FILE = "{system_name}.log"
    OUTPUT_FILE = "{system_name}.out"
    ERROR_FILE = "{system_name}.err"

    CHANGE_PERMISSION_COMMAND = "chmod 777 {path}"

    SUBMIT_COMMAND = "condor_submit {job_script_file}"
        
    STRIPPED_RECEPTOR_FILE = "stripped_receptor.pdb"
    LIGAND_FILE = "autobuilding_ligand.cif"
    EVENT_MTZ_FILE = "event.mtz"
    GRAFTED_MTZ_FILE = "grafted.mtz"
    RHOFIT_DIR = "rhofit"
    RHOFIT_EVENT_DIR = "rhofit_{}"
    RHOFIT_NORMAL_DIR = "rhofit_normal"
    RHOFIT_RESULTS_FILE = "results.txt"
    RHOFIT_RESULT_JSON_FILE = "result.json"
    RHOFIT_NORMAL_RESULT_JSON_FILE = "result_normal.json"
    RHOFIT_BEST_MODEL_FILE = "best.pdb"
    RSCC_TABLE_FILE = "rscc_table.csv"

    BUILD_DIR_PATTERN = "{pandda_name}_{dtag}_{event_idx}"


    CUMULATIVE_HITS_PLOT_FILE = "cumulative_hits.png"

    PANDDA_ANALYSES_DIR = "analyses"
    PANDDA_ANALYSE_EVENTS_FILE = "pandda_analyse_events.csv"
    PANDDA_ANALYSE_SITES_FILE = "pandda_analyse_sites.csv"
    PANDDA_PROCESSED_DATASETS_DIR = "processed_datasets"
    PANDDA_MODELLED_STRUCTURES_DIR = "modelled_structures"
    PANDDA_LIGAND_FILES_DIR = "ligand_files"
    PANDDA_PDB_FILE = "{}-pandda-input.pdb"
    PANDDA_MTZ_FILE = "{}-pandda-input.mtz"

    PANDDA_LIGAND_CIF_FILE = "ligand.cif"
    PANDDA_LIGAND_PDB_FILE = "ligand.pdb"
    PANDDA_LIGAND_SMILES_FILE = "ligand.smiles"

    PANDDA_INSPECT_EVENTS_PATH = "pandda_inspect_events.csv"
    PANDDA_EVENT_MAP_FILE = "{}-event_{}_1-BDC_{}_map.native.ccp4"
    PANDDA_EVENT_MODEL = "{}-pandda-model.pdb"

    PANDDA_PROTEIN_MASK_FILE = "protein_mask.ccp4"
    PANDDA_SYMMETRY_MASK_FILE = "symmetry_mask.ccp4"
    PANDDA_MEAN_MAP_FILE = "mean_{number}_{res}.ccp4"
    PANDDA_SIGMA_S_M_FILE = "sigma_s_m_{number}_{res}.ccp4"

    PANDDA_Z_MAP_FILE = "{dtag}-z_map.native.ccp4"
    PANDDA_EVENT_MAP_FILE = "{dtag}-event_{event_idx}_1-BDC_{bdc}_map.native.ccp4"

    PANDDA_LOG_FILE = "pandda_log.json"
    
    ELBOW_COMMAND = "cd {event_dir.event_dir}; phenix. {smiles_file.smiles_file} --output=\"{autobuilding_ligand}\""
    FINAL_LIGAND_CIF_FILE = "final_ligand_cif.cif"
    
    MASKED_PDB_FILE = "masked.pdb"
    
    RHOFIT_LOG_FILE = "rhofit.log"
    RHOFIT_OUTPUT_FILE = "rhofit.out"
    RHOFIT_ERROR_FILE = "rhofit.err"
    RHOFIT_JOB_SCRIPT_FILE = "rhofit.job"
    
    RHOFIT_SCRIPT = (
        "#!/bin/bash \n"
        "conda activate pandda \n"
        "source ~/.bashrc \n"
        "source /data/share-2/conor/xtal_software/ccp4-7.1/bin/ccp4.setup-sh \n"
        "source /data/share-2/conor/xtal_software/phenix/phenix-1.18.2-3874/phenix_env.sh \n"
        "source /data/share-2/conor/xtal_software/buster-2.10/setup.sh \n"
        "rhofit -m {mtz} -p {pdb} -l {ligand} -d {out_dir_path} -allclusters -use_2fofc -thorough "
    )
    RHOFIT_SCRIPT_FILE = "run_rhofit.sh"
    
    PHASE_GRAFTED_MTZ_FILE = "phase_grafted_mtz.mtz"
    
    RHOFIT_HIT_REGEX = "(Hit_[^\s]+)[\s]+[^\s]+[\s]+[^\s]+[\s]+([^\s]+)"
    RHOFIT_CLUSTER_BUILD_REGEX = "Hit_([^_]+)_[^_]+_([^_]+).pdb"
    
    RHOFIT_EVENT_BUILD_RESULT_FILE = "result.json"
    
    CLUSTER = "CONDOR"
    CONDOR_JOB_ID_REGEX = r"[0-9]+\."
    CONDOR_STATUS_COMMAND = r"condor_q"