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
    build_file: Path
    event_map_file: Path
    output_dir: Path
    debug: int
    
    @staticmethod
    def from_args(args: Any):
        build_file = Path(args.reference_structure_file)
        event_map_file = Path(args.structure_file)
        output_dir = Path(args.structure_file)
        debug: int = int(args.debug)
        
        return Args(
            build_file,
            event_map_file,
            output_dir,
            debug,
        )
    
    @staticmethod
    def from_cmd():
        parser = argparse.ArgumentParser()
        parser.add_argument("--build_file",
                            )
        parser.add_argument("--event_map_file",
                            )
        parser.add_argument("--output_dir",
                            )
        parser.add_argument("--debug",
                            default=2,
                            )
        args = parser.parse_args()
        
        return Args.from_args(args)

if __name__ == "__main__":
    args: Args = Args.from_cmd()
    
    # Summarise
    print(
        (
            f"Loading build from: {args.build_file}\n"
            f"Loading xmap from {args.event_map_file}\n"
        )
        )
    
    
    # Load xmap
    xmap = Xmap.from_file(args.event_map_file)
    
    # Load build
    build: Structure = Structure.from_model_path(args.build_file)
    
    # Get skeleton
    skeleton_score: xlib.SkeletonScore = xlib.SkeletonScore.from_build(
        build,
        xmap,
        contour=1.0,
    )
    
    # Print summary
    print(f"Skeleton score is: {skeleton_score}")
    
    # Save skeleton score
    json_file = args.output_dir / Constants.SKELETON_STRUCTURE_JSON_FILE
    print(f"Saving result to: {json_file}")
    
    record_dict = {xlib.Constants.SKELETON_XMAP_FILE_RECORD: str(args.event_map_file),
     xlib.Constants.SKELETON_STRUCTURE_FILE_RECORD: str(args.build_file),
     xlib.Constants.SKELETON_SCORE_RECORD: skeleton_score.skeleton_score,
     }
    with open(str(json_file), "w") as f:
        json.dump(record_dict, f)
