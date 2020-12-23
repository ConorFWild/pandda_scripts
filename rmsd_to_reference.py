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
    reference_structure_file: Path
    structure_file: Path
    debug: int
    
    @staticmethod
    def from_args(args: Any):
        reference_structure_file = Path(args.reference_structure_file)
        structure_file = Path(args.structure_file)
        debug: int = int(args.debug)
        
        return Args(
            reference_structure_file,
            structure_file,
            debug,
        )
    
    @staticmethod
    def from_cmd():
        parser = argparse.ArgumentParser()
        parser.add_argument("--reference_structure_file",
                            )
        parser.add_argument("--structure_file",
                            )
        parser.add_argument("--debug",
                            default=2,
                            )
        args = parser.parse_args()
        
        return Args.from_args(args)

if __name__ == "__main__":
    args: Args = Args.from_cmd()
    
    reference_structure = Structure.from_model_path()
    
    build_structure = Structure.from_model_path()
    
    rmsd: RMSD = RMSD.from_structures(
        reference_structure=reference_structure,
        build_structure=build_structure,
    )