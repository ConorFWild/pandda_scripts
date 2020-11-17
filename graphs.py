from __future__ import annotations
from batch_pandda import OUTPUT_FILE

from os import system, write, mkdir
import shutil
from sys import path, stderr, stdin
from typing import *
import dataclasses

from pathlib import Path
import argparse
import subprocess

import numpy as np
import pandas as pd

import joblib

import gemmi

import seaborn as sns

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
        "source ~/.bashrc \n"
        "conda activate pandda \n"
        "rhofit -m {mtz} -p {pdb} -l {ligand} -d {out_dir_path} -allclusters -use_2fofc -thorough"
    )
    RHOFIT_SCRIPT_FILE = "run_rhofit.sh"
    
    PHASE_GRAFTED_MTZ_FILE = "phase_grafted_mtz.mtz"
    
    DATASET_PDB_FILE = "final.pdb"
    DATASET_MTZ_FILE = "final.mtz"
    
    DATASET_RESOLUTION_GRAPH = "dataset_resolution.png"
    EVENT_SIZE_GRAPH = "event_size.png"



@dataclasses.dataclass()
class Args:
    data_dirs_dir: Path
    pandda_dirs_dir: Path
    autobuild_dirs_dir: Path
    graph_dir: Path
    debug: int
    
    @staticmethod
    def from_args(args: Any):
        data_dirs_dir = Path(args.data_dirs_dir)
        pandda_dirs_dir = Path(args.pandda_dirs_dir)
        autobuild_dirs_dir = Path(args.autobuild_dirs_dir)
        graph_dir = Path(args.graph_dir)
        debug: int = int(args.debug)
        
        return Args(
            data_dirs_dir,
            pandda_dirs_dir,
            autobuild_dirs_dir,
            graph_dir,
            debug,
        )
    
    @staticmethod
    def from_cmd():
        parser = argparse.ArgumentParser()
        parser.add_argument("--data_dirs_dir",
                            )
        parser.add_argument("--pandda_dirs_dir",
                            )
        parser.add_argument("--autobuild_dirs_dir",
                            )
        parser.add_argument("--graph_dir",
                            default=2,
                            )
        parser.add_argument("--debug",
                            default=2,
                            )
        args = parser.parse_args()
        
        return Args.from_args(args)
    


@dataclasses.dataclass()
class System:
    system: str
    
    def __hash__(self) -> int:
        return hash(self.system)

@dataclasses.dataclass()
class Dtag:
    dtag: str
    
    def __hash__(self) -> int:
        return hash(self.dtag)

@dataclasses.dataclass()
class EventIDX:
    event_idx: str
    
    def __hash__(self) -> int:
        return hash(self.event_idx)

@dataclasses.dataclass()
class EventID:
    system: System
    dtag: Dtag
    event_idx: EventIDX
    
    def __hash__(self) -> int:
        return hash(
            (self.system,
                     self.dtag,
                     self.event_idx,
                     )
                    )
    
@dataclasses.dataclass()
class Event:
    event_input_dir: Path
    event_output_dir: Path
    system: System
    dtag: Dtag
    event_idx: EventIDX
    bdc: float
    x: float
    y: float
    z: float
    size: int
    
@dataclasses.dataclass()
class DatasetID:
    system: System
    dtag: Dtag
    
    def __hash__(self) -> int:
        return hash(
            (self.system,
                     self.dtag,
                     )
                    )
    
@dataclasses.dataclass()
class Dataset:
    dataset_id: DatasetID
    resolution: float    

# ##############
# Functions
# ##############

def get_event_input_dir(system: System, dtag: Dtag, event_idx: EventIDX, pandda_dirs_dir: Path) -> Path:
    return pandda_dirs_dir / system.system /Constants.PANDDA_PROCESSED_DATASETS_DIR / dtag.dtag

def get_event_output_dir(system: System, dtag: Dtag, event_idx: EventIDX, autobuild_dirs_dir: Path) -> Path:
    return autobuild_dirs_dir / system.system / dtag.dtag / event_idx.event_idx
    
def get_event_table_dict(path_list: List[Path]) -> Dict[System, pd.DataFrame]:
    event_table_dict = {}
    
    for path in path_list:    
        system: System = System(path.name)
        
        event_table_file: Path = path / Constants.PANDDA_ANALYSES_DIR / Constants.PANDDA_ANALYSE_EVENTS_FILE
        
        if Constants.DEBUG >0: print(event_table_file)
        
        if event_table_file.exists():
            # try:
            event_table: pd.DataFrame = pd.read_csv(str(event_table_file))
            event_table_dict[system] = event_table
            # except Exception as e:
            #     print(e)
            #     print(f"event_table_file seems empty?")
            #     continue
    
    return event_table_dict

def get_event_id(system: System, row: pd.Series) -> EventID:
    dtag: Dtag = Dtag(row["dtag"])
    event_idx: EventIDX = EventIDX(row["event_idx"])
    return EventID(
        system=system,
        dtag=dtag,
        event_idx=event_idx,
        )
    
def get_event(system: System, row: pd.Series, pandda_dirs_dir: Path, autobuild_dirs_dir: Path) -> Event:

    dtag: Dtag = Dtag(row["dtag"])
    event_idx: EventIDX = EventIDX(str(row["event_idx"]))
    bdc = row["1-BDC"]
    x = row["x"]
    y = row["y"]
    z = row["z"]
    
    
    event_input_dir: Path = get_event_input_dir(system, dtag, event_idx, pandda_dirs_dir)
    event_output_dir: Path = get_event_output_dir(system, dtag, event_idx, autobuild_dirs_dir)
    
    size: int = int(row["cluster_size"])
    
    return Event(
        event_input_dir,
        event_output_dir,
        system=system,
        dtag=dtag,
        event_idx=event_idx,
        bdc=bdc,
        x=x,
        y=y,
        z=z,
        size=size,
    )
    

def get_event_dict(event_table_dict: Dict[System, pd.DataFrame], pandda_dirs_dir: Path, autobuild_dirs_dir: Path) -> Dict[EventID, Event]:
    event_dict: Dict[EventID, Event] = {}
    
    for system, event_table in event_table_dict.items():
        for index, row in event_table.iterrows():
            event_id: EventID = get_event_id(system, row)
            event: Event = get_event(system, row, pandda_dirs_dir, autobuild_dirs_dir)
            
            event_dict[event_id] = event
            
    return event_dict

def get_dataset_dict(system_path_list: List[Path]) -> Dict[DatasetID, Dataset]:
    dataset_dict: Dict[DatasetID, Dataset] = {}
    
    for system_path in system_path_list:
        system: System = System(system_path.name)
        
        dataset_dir_list: List[Path] = list(system_path.glob("*"))
        
        for dataset_dir in dataset_dir_list:
            dtag: Dtag = Dtag(dataset_dir.name)
            dataset_pdb_file: Path = dataset_dir / Constants.DATASET_PDB_FILE
            dataset_mtz_file: Path = dataset_dir / Constants.DATASET_MTZ_FILE
            
            dataset_mtz: gemmi.Mtz = gemmi.read_mtz_file(str(dataset_mtz_file))
            dataset_resolution: float = dataset_mtz.resolution_high()
            
            dataset_id: DatasetID = DatasetID(system, dtag)
            
            dataset: Dataset = Dataset(dataset_id, dataset_resolution)
            
            dataset_dict[dataset_id] = dataset
    
    return dataset_dict

# #########
# Graphs
# ###########
def make_resolution_graph(
        dataset_dict: Dict[DatasetID, Dataset],
        dataset_resolution_graph_file: Path,
        ):
    resolution_list: List[float] = [dataset.resolution for dataset in dataset_dict.values()]
    plot = sns.displot(resolution_list)
    fig = plot.get_figure()
    fig.savefig(str(dataset_resolution_graph_file))
    

# ###########
# Script
# ###########


def main():
    # get args
    args: Args = Args.from_cmd()
    if args.debug > 0:
        Constants.DEBUG = args.debug

    # Get system directories
    system_path_list: List[Path] = list(args.data_dirs_dir.glob("*"))
    
    # Get datasets
    dataset_dict: Dict[DatasetID, Dataset] = get_dataset_dict(system_path_list)

    # Get pandda directories
    system_path_list: List[Path] = list(path for path in args.pandda_dirs_dir.glob("*") if path.is_dir())

    # Get events tables: List[Path] -> Dict[EventID, Event]
    event_table_dict: Dict[System, pd.DataFrame] = get_event_table_dict(system_path_list)
    
    # Get events
    event_dict: Dict[EventID, Event] = get_event_dict(event_table_dict, args.pandda_dirs_dir, args.autobuild_dirs_dir)

    # Print num datasets
    print(f"Number of datasets is: {len(dataset_dict)}")

    # Make resolution graph
    make_resolution_graph(
        dataset_dict,
        args.graph_dir / Constants.DATASET_RESOLUTION_GRAPH,
        )
    
    # Make event size
    make_event_size_graph(
        event_dict,
        args.graph_dir / Constants.EVENT_SIZE_GRAPH,
    )    
        

if __name__ == "__main__":
    main()