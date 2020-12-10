from __future__ import annotations
from batch_rhofit import get_event_table_dict
from re import match
from batch_pandda import OUTPUT_FILE

from os import system, write, mkdir
import shutil
from sys import path, stderr, stdin, stdout
from typing import *
import dataclasses

from pathlib import Path
import argparse
import subprocess
import re
import json
import time

import numpy as np
import pandas as pd

import joblib

import seaborn as sns

import gemmi

from xlib import *

def tail(file: Path, n: int = 20) -> str:
    command: str = f"tail -n {n} {str(file)}"
    print(command)
    p = subprocess.Popen(command,
                         shell=True,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE,
                         )
    stdout, stderr = p.communicate()
    return str(stdout)

def errored(string: str) -> bool:
    pattern: Pattern = re.compile("Traceback")
    matches: List[Any] = re.findall(pattern, string)
    if len(matches) > 0:
        return True
    else:
        return False

def make_event_distance_graph(event_distance_dict: Dict[Dtag, float], path: Path):
    
    def categorise(number: float) -> str:
        if number< 0:
            return "Missed!"
        if number < 2:
            return "<2"
        elif number <5:
            return "2-5"
        else:
            return ">5"
    
    data = {
        "Dtag": [dtag for dtag in event_distance_dict],
        "Distance": [categorise(event_distance_dict[dtag]) for dtag in event_distance_dict],
    }
    
    table = pd.DataFrame.from_dict(data)

    g = sns.catplot(x="Distance",
                    kind="count",
                    data=table,
                    )

    
    g.savefig(str(path))
    
def make_successful_pandda_plot(
    system_path_dict, 
    event_table_dict, 
    err_dict, 
    path: Path):
    
    def categorise(system, event_table_dict, err_dict) -> str:
        if system in event_table_dict:
            return "Complete"
        elif system in err_dict:
            return "Errored"
        else:
            return "Unknown"
        
    data = {
        "System": [system.system for system in system_path_dict],
        "Status": [categorise(system, event_table_dict, err_dict) 
                   for system
                   in system_path_dict
                   ],
    }
    table = pd.DataFrame.from_dict(data)

    g = sns.catplot(x="Status",
                    kind="count",
                    data=table,
                    )
    
    g.savefig(str(path))
    
def make_event_plot(
    system_path_dict, 
    event_table_dict, 
    reference_structure_dict: ReferenceStructureDict,
    event_dict: EventDict,
    path: Path):
    
    def categorise(system, event_table_dict, err_dict) -> str:
        if system in event_table_dict:
            return "Complete"
        elif system in err_dict:
            return "Errored"
        else:
            return "Unknown"
        
    def get_num_events(system: System, event_dict: EventDict) -> int:
        system_events = [event_id for event_id in event_dict if event_id.system.system == system.system]
        return len(system_events)
        
    def get_num_datasets(system: System, system_path_dict: SystemPathDict) -> int:
        processed_dataset_path: Path = system_path_dict[system] / Constants.PANDDA_PROCESSED_DATASETS_DIR
        processed_dataset_dir_list: List[Path] = list(processed_dataset_path.glob("*"))
        return len(processed_dataset_dir_list)
        
    def get_num_hits(system: System, reference_structure_dict: ReferenceStructureDict,
                     event_dict: EventDict) -> int:
        
        dtag_structures = []
        
        for dtag in reference_structure_dict:
            dtag_system = System.from_dtag(dtag.dtag)
            if not dtag_system.system == system.system:
                continue
            
            # Check if an event
            dtag_events = [event_id for event_id in event_dict if event_id.dtag.dtag == dtag.dtag]
            
            if len(dtag_events) > 0:
                dtag_structures.append(dtag)
                
        return len(dtag_structures)
    
    def get_num_target_hits(system: System, reference_structure_dict: ReferenceStructureDict) -> int:
        dtag_structures = []
        
        for dtag in reference_structure_dict:
            dtag_system = System.from_dtag(dtag.dtag)
            if not dtag_system.system == system.system:
                continue
            dtag_structures.append(dtag)
        
        return len(dtag_structures)
            
            
    data = []    
    # Make events/datasets
    for system in system_path_dict:
        num_events: int = get_num_events(system, event_table_dict) 
        num_datasets: int = get_num_datasets(system, system_path_dict)
        
        # Get percentage
        if num_datasets == 0:
            percentage = 0.0
        else:
            percentage = num_events / num_datasets
            
        record = {
            "System": system.system,
            "Kind": "num events/num datasets",
            "Percentage": percentage,
        }
        data.append(record)
        
    # Make hits/events
    for system in system_path_dict:
        num_hits: int =  get_num_hits(system, reference_structure_dict, event_dict,)
        num_events: int = get_num_events(system, event_table_dict) 
        
        print(num_hits)
        print(num_events)
        
        # Get percentage
        if num_events == 0:
            percentage = 0.0
        else:
            percentage = num_hits / num_events
        
        record = {
            "System": system.system,
            "Kind": "captured hits/number of events",
            "Percentage": percentage,
        }
        data.append(record)
    
    # makes hits / dataset hits
    for system in system_path_dict:
        num_hits:int = get_num_hits(system, reference_structure_dict, event_dict) 
        num_target_hits:int = get_num_target_hits(system, reference_structure_dict)
        
        # Get percentage
        if num_target_hits == 0:
            percentage = 0.0
        else:
            percentage = num_hits / num_target_hits
        
        record = {
            "System": system.system,
            "Kind": "captured hits/target hits",
            "Percentage": percentage,
        }
        data.append(record)
    
    
    table = pd.DataFrame(data)

    g = sns.catplot(x="System",
                    y="Percentage",
                    hue="Kind",
                    kind="bar",
                    data=table,
                    aspect=12,
                    )
    g.set_xticklabels(rotation=30)
    
    g.savefig(str(path))

@dataclasses.dataclass()
class Args:
    pandda_dirs_dir: Path
    autobuild_dirs_dir: Path
    graph_dir: Path
    debug: int
    
    @staticmethod
    def from_args(args: Any):
        pandda_dirs_dir = Path(args.pandda_dirs_dir)
        autobuild_dirs_dir = Path(args.autobuild_dirs_dir)
        graph_dir: Path = Path(args.graph_dir)
        debug: int = int(args.debug)
        
        return Args(
            pandda_dirs_dir,
            autobuild_dirs_dir,
            graph_dir,
            debug,
        )
    
    @staticmethod
    def from_cmd():
        parser = argparse.ArgumentParser()
        parser.add_argument("--pandda_dirs_dir",
                            )
        parser.add_argument("--autobuild_dirs_dir",
                            )
        parser.add_argument("--graph_dir",
                            )
        parser.add_argument("--debug",
                            default=2,
                            )
        args = parser.parse_args()
        
        return Args.from_args(args)

def main():
    # Get args
    args: Args = Args.from_cmd()
    
    # Get systems
    pandda_system_path_dict: SystemPathDict = SystemPathDict.from_dir(args.pandda_dirs_dir)
    if args.debug > 0: print(f"Found {len(pandda_system_path_dict)} systems")

    # Get PanDDAs
    event_table_dict: EventTableDict = EventTableDict.from_system_path_dict(pandda_system_path_dict)
    if args.debug > 0: print(f"Found {len(event_table_dict)} system tables")    
    
    # Get events
    event_dict: EventDict = EventDict.from_event_tables(event_table_dict, args.pandda_dirs_dir, args.autobuild_dirs_dir)
    if args.debug > 0: 
        system_dict: Dict[System, None] = {build_id.system: None for build_id in event_dict}
        print(f"Got events for {len(system_dict)} systems")
        dtag_dict: Dict[Dtag, None] = {build_id.dtag: None for build_id in event_dict}
        print(f"Got events for {len(dtag_dict)} datasets")
        print(f"Found {len(event_dict)} events")
    
    # Unfinished systems dict
    unfinished_system_path_dict: SystemPathDict = SystemPathDict(
        {system: pandda_system_path_dict[system]
         for system
         in pandda_system_path_dict
         if not (pandda_system_path_dict[system] /Constants.PANDDA_ANALYSES_DIR /Constants.PANDDA_ANALYSE_EVENTS_FILE).exists()
         }
        )
    if args.debug > 0: print(f"Found {len(unfinished_system_path_dict)} unfinished systems")    

    
    # err dict
    err_dict: Dict[System, Path] = {
        system: args.pandda_dirs_dir / Constants.PANDDA_JOB_ERROR_FILE.format(system_name=system.system)
        for system
        in unfinished_system_path_dict
        }
    if args.debug > 0: 
        for system in err_dict:
            print((
                f"# # {system.system} \n"
                f"{err_dict[system]} \n"
                f"{tail(err_dict[system])} \n"
            )
                  )    
        # Print num errored systems
        errored_target_list: List[System] = [system 
                                                    for system 
                                                    in err_dict 
                                                    if errored(tail(err_dict[system]))
                                                    ]
        print(f"Number of systems errored: {len(errored_target_list)}")
        print(f"Errored systems: {errored_target_list}")

    # Pandda status graph
    make_successful_pandda_plot(pandda_system_path_dict,
                                event_table_dict,
                                err_dict,
                                args.graph_dir / "pandda_gemmi_status_graph.png",
                                )

    # Get references
    reference_structure_dict: ReferenceStructureDict = ReferenceStructureDict.from_system_path_dict(pandda_system_path_dict)
    if args.debug > 0: 
        print(f"Found {len(reference_structure_dict)} reference structures")    

    # events
    if args.debug > 0: 
        make_event_plot(
            pandda_system_path_dict, 
            event_table_dict, 
            reference_structure_dict,
            event_dict,
            args.graph_dir / "pandda_gemmi_all_event_summary.png",
            )
        
    # load jsons
    # for system_path in system_paths:
    #     json_path = system_path / "log.json"
    #     with open(str(json_path)):
    #         json_dict = json.load(str(json_path))

    # Check for event tables
    event_distance_dict: Dict[Dtag, float] = get_event_distance_from_reference_model_dict(
        event_dict,
        reference_structure_dict,
        pandda_system_path_dict,
    )
    if args.debug > 0: 
        print(event_distance_dict)
        print(f"Found {len(event_distance_dict)} closest events to know hits")    
    make_event_distance_graph(event_distance_dict, 
                              args.graph_dir / "pandda_gemmi_all_closest_event_distance.png",
                              )

    # events by 
    
    # Output

if __name__ == "__main__":
    main()