from typing import *
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

JOB_SCRIPT_FILE = "{system_name}.job"

LOG_FILE = "{system_name}.log"
OUTPUT_FILE = "{system_name}.out"
ERROR_FILE = "{system_name}.err"


PANDDA_SCRIPT = """
#!/bin/bash

source ~/.bashrc
conda activate pandda

pandda=\""/data/share-2/conor/pandda/code/pandda_gemmi/pandda_gemmi/analyse.py"\"
data_dirs=\"{data_dirs}\"
cluster_cutoff_distance_multiplier=1.1 
out_dir=\"{out_dir}\"
pdb_regex=\"final.pdb\"
mtz_regex=\"final.mtz\"

python $pandda --data_dirs=$data_dirs --cluster_cutoff_distance_multiplier=$cluster_cutoff_distance_multiplier --min_blob_volume=$min_blob_volume --out_dir=$out_dir --pdb_regex=$pdb_regex --mtz_regex=$mtz_regex


"""


PANDDA_SCRIPT_FILE = "{system_name}.sh"

SUBMIT_COMMAND = "condor_submit {job_file}"

@dataclasses.dataclass()
class Args:
    data_dirs_dir: Path
    out_dirs_dir: Path
    request_memory: int
    debug: bool
    
    @staticmethod
    def from_args(args: Any):
        data_dirs_dir = Path(args.data_dirs_dir)
        out_dirs_dir = Path(args.out_dirs_dir)
        request_memory= int(args.request_memory)
        debug: bool = bool(args.debug)
        
        return Args(
            data_dirs_dir,
            out_dirs_dir,
            request_memory,
            debug,
        )
    

def main():
    # get args
    parser = argparse.ArgumentParser()
    parser.add_argument("--data_dirs_dir",
                        )
    parser.add_argument("--out_dirs_dir",
                        )
    parser.add_argument("--request_memory",
                        )
    parser.add_argument("--debug",
                        )
    args = Args.from_args(parser.parse_args())

    # Get data dirs
    system_path_list = list(args.data_dirs_dir.glob("*"))
    
    # Make pandda commands
    pandda_script_dict = {}
    for system_path in system_path_list:
        system_name = system_path.name
        data_dirs = str(system_path)
        out_dir = args.out_dirs_dir / system_name
        pandda_script = PANDDA_SCRIPT.format(
            data_dirs=data_dirs,
            out_dir=out_dir,
        )
        pandda_script_dict[system_name] = pandda_script
        
        if args.debug:
            print(f"# # System: {system_name}: pandda script")
            print(pandda_script)
            break

    # Make command scripts
    pandda_script_files = {}
    for system_name, pandda_script in pandda_script_dict.items():
        pandda_script_file = args.out_dirs_dir / PANDDA_SCRIPT_FILE.format(system_name=system_name)
        pandda_script_files[system_name] = pandda_script_file
        
        with open(pandda_script_file, "w") as f:
            f.write(pandda_script)
            
        if args.debug:
            print(f"# # System: {system_name}: pandda script file")
            print(pandda_script_file)
            break
        
    # Make job scripts
    job_script_dict = {}
    for system_name, pandda_script_file in pandda_script_files.items():
        log_file = LOG_FILE.format(system_name=system_name)
        output_file = OUTPUT_FILE.format(system_name=system_name)
        error_file = ERROR_FILE.format(system_name=system_name)
        
        job_script = JOB_SCRIPT.format(
            pandda_script_file=pandda_script_file,
            log_file=log_file,
            output_file=output_file,
            error_file=error_file,
            request_memory=args.request_memory,
        )
        
        job_script_dict[system_name] = job_script
        
        if args.debug:
            print(f"# # System: {system_name}: Job script")
            print(job_script)
            break
    
    
    # Writer job files
    job_script_file_dict = {}
    for system_name, job_script in job_script_dict.items():
        job_script_file = args.out_dirs_dir / JOB_SCRIPT_FILE.format(system_name=system_name)
    
        with open(str(job_script_file), "w") as f:
            f.write(job_script)
    
        job_script_file_dict[system_name] = job_script_file
        
        if args.debug:
            print(f"# # System: {system_name}: job script file")
            print(job_script_file)
            break
        
            
    # Make Submit commands
    submit_command_dict = {}
    for system_name, pandda_script_file in pandda_script_files.items():
        submit_command = SUBMIT_COMMAND.format(pandda_script_file=pandda_script_file)
        submit_command_dict[system_name] = submit_command
        
        if args.debug:
            print(f"# # System: {system_name}: submit command")
            print(submit_command)
            break
        

    
    #     Submit
    for system_name, command in submit_command_dict.items():
        p = subprocess.Popen(command,
                         shell=True,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE,
                         )
        stdout, stderr = p.communicate()
        
        if args.debug:
            print(f"# # System: {system_name}: submit stdout, stderr")
            print(stdout)
            print(stderr)
            break
        
                

if __name__ == "__main__":
    main()