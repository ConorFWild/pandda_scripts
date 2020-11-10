import typing
import dataclasses

from pathlib import Path
import argparse
import subprocess
import json

def main():
    # Get args
    parser = argparse.ArgumentParser()
    parser.add_argument("--out_dirs",
                        )
    args = parser.parse_args()
    
    # Get pandda dirs
    system_paths = list(args.out_dirs.glob("*"))
    
    # load jsons
    for system_path in system_paths:
        json_path = system_path / "log.json"
        with open(str(json_path))
        json_dict = json.load(str(json_path))
    
    # check which are finished
    
    # Output

if __name__ == "__main__":
    main()