
from __future__ import annotations
from lib import EventTableDict
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
    if Constants.DEBUG > 0: print(f"Found {len(pandda_system_path_dict)} systems")

    # Get PanDDAs
    event_table_dict: EventTableDict = EventTableDict.from_system_path_dict(pandda_system_path_dict,)
    if Constants.DEBUG > 0: print(f"Found {len(event_table_dict)} system tables")
    
    # Get builds
    build_dict: BuildDict = BuildDict.from_autobuild_dir(args.autobuild_dirs_dir)
    if Constants.DEBUG > 0: print(f"Found {len(build_dict)} builds")
    
    # Get references
    reference_structure_dict: ReferenceStructureDict = ReferenceStructureDict.from_system_path_dict(pandda_system_path_dict)
    if Constants.DEBUG > 0: print(f"Found {len(reference_structure_dict)} reference structures")
    
    # Get rmsds
    rmsd_dict: RMSDDict = RMSDDict.from_build_dict(build_dict, reference_structure_dict)
    if args.debug > 0: print(rmsd_dict)
    
    # Get rsccs
    # rscc_dict: RSCCDict = RSCCDict.from_build_dict(build_dict)
    
    # 

    
    
if __name__ == "__main__":
    main()