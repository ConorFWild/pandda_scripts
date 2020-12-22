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

if __name__ == "__main__":
    args: Args = Args.from_cmd()
    