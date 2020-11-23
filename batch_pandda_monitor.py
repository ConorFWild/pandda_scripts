from __future__ import annotations
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

import gemmi

from xlib import *

def tail(file: Path, n: int = 20) -> str:
    command: str = f"tail -n {n} {str(path)}"
    p = subprocess.Popen(command,
                         shell=True,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE,
                         )
    stdout, stderr = p.communicate()
    return str(stdout)


@dataclasses.dataclass()
class Args:
    pandda_dirs_dir: Path
    autobuild_dirs_dir: Path
    debug: int
    
    @staticmethod
    def from_args(args: Any):
        pandda_dirs_dir = Path(args.pandda_dirs_dir)
        autobuild_dirs_dir = Path(args.autobuild_dirs_dir)
        debug: int = int(args.debug)
        
        return Args(
            pandda_dirs_dir,
            autobuild_dirs_dir,
            debug,
        )
    
    @staticmethod
    def from_cmd():
        parser = argparse.ArgumentParser()
        parser.add_argument("--pandda_dirs_dir",
                            )
        parser.add_argument("--autobuild_dirs_dir",
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
    
    # Unfinished systems dict
    unfinished_system_path_dict: SystemPathDict = SystemPathDict(
        {system: pandda_system_path_dict[system]
         for system
         in pandda_system_path_dict
         if not (pandda_system_path_dict[system] /Constants.PANDDA_ANALYSES_DIR /Constants.PANDDA_ANALYSE_EVENTS_FILE).exists()
         }
        )
    
    # err dict
    err_dict: Dict[System, Path] = {
        system: args.pandda_dirs_dir / Constants.PANDDA_JOB_ERROR_FILE.format(system.system)
        for system
        in unfinished_system_path_dict
        }
    if args.debug > 0: 
        for system in err_dict:
            print((
                f"# # {system.system} \n"
                f"{tail(err_dict[system])} \n"
            )
                  )    

    
    # out dict
    
    # load jsons
    # for system_path in system_paths:
    #     json_path = system_path / "log.json"
    #     with open(str(json_path)):
    #         json_dict = json.load(str(json_path))

    # Check for event tables
    
    # check which are finished
    
    # Output

if __name__ == "__main__":
    main()