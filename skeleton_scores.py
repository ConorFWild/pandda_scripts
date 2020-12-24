from __future__ import annotations
from re import match

import os
import shutil
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
    out_dir: Path
    
    @staticmethod
    def from_cmd():
        
        fields: Tuple[dataclasses.Field, ...] = dataclasses.fields(Args)
        parser = argparse.ArgumentParser()
        for field in fields:
            parser.add_argument(f"--{field.name}")
        args = parser.parse_args()
        
        args_dict = vars(args)
        
        typed_args = [field.type(args_dict[field.name]) for field in fields]
        
        return Args.from_args(args)

def make_output_dirs(out_dir: Path, build_dict: BuildDict):
    
    for build_id in build_dict:
        
        # System
        system_dir = out_dir / build_id.system.system
        if not system_dir.exists():
            os.mkdir(str(system_dir))
        
        # Dtag
        dtag_dir = system_dir / build_id.dtag.dtag
        if not dtag_dir.exists():
            os.mkdir(str(dtag_dir))
        
        # event_idx
        event_dir = dtag_dir / build_id.event_idx.event_idx
        if not event_dir.exists():
            os.mkdir(str(event_dir))
        
        # build cluster
        cluster_dir = event_dir / str(build_id.build_cluster.build_cluster_id)
        if not cluster_dir.exists():
            os.mkdir(str(cluster_dir))
        
        # build number
        build_dir = cluster_dir / str(build_id.build_number.build_number_id)
        if not build_dir.exists():
            os.mkdir(str(build_dir))
        
    

    
if __name__ == "__main__":
    args: Args = Args.from_cmd()
    
    # Summarise options
    print(
        (
            f"pandda_dirs_dir: {args.pandda_dirs_dir}"
            f"Autobuilds dir: {args.autobuild_dirs_dir}"
            f"Output dir: {args.out_dir}"
        )
    )

    pandda_system_path_dict: SystemPathDict = SystemPathDict.from_dir(args.pandda_dirs_dir)
    if args.debug > 0: print(f"Found {len(pandda_system_path_dict)} systems")

    # Get PanDDAs
    event_table_dict: EventTableDict = EventTableDict.from_system_path_dict(pandda_system_path_dict,)
    if args.debug > 0: print(f"Found {len(event_table_dict)} system tables")
    
    # Get references
    reference_structure_dict: ReferenceStructureDict = ReferenceStructureDict.from_system_path_dict(pandda_system_path_dict)
    if args.debug > 0: 
        print(f"Found {len(reference_structure_dict)} reference structures")    
    
    # Get events
    event_dict: EventDict = EventDict.from_event_tables(event_table_dict, args.pandda_dirs_dir, args.autobuild_dirs_dir)
    
    # Get builds
    build_dict: BuildDict = BuildDict.from_autobuild_dir(args.out_dir, args.autobuild_dirs_dir)
    
    # Make output dirs
    make_output_dirs(build_dict)
    
    # Run Skeleton scoring
    for build_id in build_dict:
        # Get out dir
        out_dir: Path = args.out_dir / build_id.system.system / build_id.dtag.dtag / build_id.event_idx.event_idx / str(build_id.build_cluster.build_cluster_id) / str(build_id.build_number.build_number_id)
        
        # Get build file
        build: Build = build_dict[build_id]
        build_file: str = str(build.build_file)
        
        # Get event map file
        event_id = EventID(build_id.system, build_id.dtag, build_id.event_idx)
        event: Event = event_dict[event_id]
        # event_dir: Path = pandda_system_path_dict[event_id.system] / xlib.Constants.PANDDA_PROCESSED_DATASETS_DIR / build_id.dtag.dtag
        event_dir: Path = event.event_input_dir
        event_map_file: Path = event_dir / xlib.Constants.PANDDA_EVENT_MAP_FILE.format(
            dtag=event.dtag.dtag,
            event_idx=event.event_idx.event_idx,
            bdc=event.bdc,
            )
        
        # Run
        env = "conda activate env_rdkit"
        command = "{env}; python --build_file {build_file} --event_map_file {event_map_file} --output_dir {out_dir}".format(
            env=env,
            build_file=build_file,
            event_map_file=event_map_file,
            out_dir=out_dir,
        )
        p = subprocess.Popen(command,
                         shell=True,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE,
                         )
        stdout, stderr = p.communicate()
        print(stdout)
        print(stderr)
        
        