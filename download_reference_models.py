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



@dataclasses.dataclass()
class Args:
    data_dirs_dir: Path
    reference_structure_dir: Path
    debug: int
    
    @staticmethod
    def from_args(args: Any):
        data_dirs_dir = Path(args.pandda_dirs_dir)
        reference_structure_dir = Path(args.autobuild_dirs_dir)
        debug: int = int(args.debug)
        
        return Args(
            data_dirs_dir,
            reference_structure_dir,
            debug,
        )
    
    @staticmethod
    def from_cmd():
        parser = argparse.ArgumentParser()
        parser.add_argument("--data_dirs_dir",
                            )
        parser.add_argument("--reference_structure_dir",
                            )
        parser.add_argument("--debug",
                            default=2,
                            )
        args = parser.parse_args()
        
        return Args.from_args(args)

if __name__ == "__main__":
    args: Args = Args.from_cmd()
    
    system_path_dict: SystemPathDict = SystemPathDict.from_dir(args.pandda_dirs_dir)
    
    reference_structure_dict: ReferenceStructureDict = ReferenceStructureDict.from_system_path_dict(system_path_dict)
    
    for dtag in reference_structure_dict:
        structure: Structure = reference_structure_dict[dtag]
        
        path: Path = args.reference_structure_dir / xlib.Constants.REFERENCE_STRUCTURE_FILE.format(dtag=dtag.dtag)
        
        structure.to_pdb(path)