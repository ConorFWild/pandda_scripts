import typing
import dataclasses

from pathlib import Path
import argparse
import subprocess

SUBMIT_COMMAND = "qsub -pe smp 10 -l m_mem_free=15G,h_vmem=200G {pandda_script_file}"

JOB_SCRIPT = """
####################                                                    
# 
# Example 1                                                            
# Simple HTCondor submit description file                                    
#                                                                       
####################    

Executable   = {pandda_script_file}                                                    
Log          = {log_file}    
Output = {output_file}
Error = {error_file}
                                                                        
request_memory = {request_memory} GB


Queue
"""

PANDDA_SCRIPT = """
#!/bin/bash

source ~/.bashrc
conda activate pandda

pandda=\""/data/share-2/conor/pandda/code/pandda_gemmi/pandda_gemmi/analyse.py"\"
data_dirs=\"{data_dirs\"
cluster_cutoff_distance_multiplier=1.1 
out_dir=\"{}\"
pdb_regex=\"final.pdb\"
mtz_regex=\"final.mtz\"


python $pandda --data_dirs=$data_dirs --cluster_cutoff_distance_multiplier=$cluster_cutoff_distance_multiplier --min_blob_volume=$min_blob_volume --out_dir=$out_dir --pdb_regex=$pdb_regex --mtz_regex=$mtz_regex


"""

PANDDA_SCRIPT_FILE = "{system_name}.sh"

@dataclasses.dataclass()
class Args:
    data_dirs_dir: Path
    out_dirs_dir: Path
    
    @staticmethod
    def from_args(args: Any):
        data_dirs_dir = Path(args.data_dirs_dir)
        out_dirs_dir = Path(args.out_dirs_dir)
        
        return Args(
            data_dirs_dir,
            out_dirs_dir
        )
    

def main():
    # get args
    parser = argparse.ArgumentParser()
    parser.add_argument("--data_dirs_dir",
                        )
    parser.add_argument("--out_dirs_dir",
                        )
    args = Args.from_args(parser.parse_args())

    # Get data dirs
    system_path_list = list(args.data_dirs_dir.glob("*"))
    
    # Make pandda commands
    commands = {}
    for system_path in system_path_list:
        system_name = system_path.name
        data_dirs = str(system_path)
        out_dir = args.out_dirs_dir / system_name
        command = PANDDA_SCRIPT.format(
            data_dirs=data_dirs,
            out_dir=out_dir,
        )
        commands[system_name] = command

    # Make command scripts
    pandda_script_files = {}
    for system_name, command in commands.items():
        pandda_script_file = args.out_dirs_dir / PANDDA_SCRIPT_FILE.format(system_name=system_name)
        pandda_script_files[system_name] = pandda_script_file
        
        with open(pandda_script_file, "w") as f:
            f.write(command)
            
    # Make Submit commands
    submit_command_dict = {}
    for system_name, pandda_script_file in pandda_script_files.items():
        submit_command = SUBMIT_COMMAND.format(pandda_script_file=pandda_script_file)
        submit_command_dict[system_name] = submit_command
    
    #     Submit
    for system_name, command in submit_command_dict.items():
        p = subprocess.Popen(command,
                         shell=True,
                         )
        p.communicate()
                

if __name__ == "__main__":
    main()