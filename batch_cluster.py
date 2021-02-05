import fire

from pathlib import Path

import os
import subprocess

class Constants:
    SCRIPT = (
        "#!/bin/bash\n"
        ". /data/share-2/conor/anaconda3/etc/profile.d/conda.sh; conda activate env_rdkit; source /data/share-2/conor/xtal_software/ccp4-7.1/bin/ccp4.setup-sh; source /data/share-2/conor/xtal_software/phenix/phenix-1.18.2-3874/phenix_env.sh; source /data/share-2/conor/xtal_software/buster-2.10/setup.sh\n"
        "python /data/share-2/conor/pandda/pandda_scripts/local_cluster_mask.py {data_dir}/{system} {out_dir}/{system} \"*.pdb\" \"*.mtz\" \"FWT,PHWT\" \"grid\"\n"
    )
    
    JOB = (
        "####################\n"
        "#\n"
        "# Example 1\n"
        "# Simple HTCondor submit description file\n"
        "#\n"
        "####################\n"
        
        "Executable   = {out_dir}/{system}.sh\n"
        "Log          = {out_dir}/{system}.log\n" 
        "Output = {out_dir}/{system}.out\n"
        "Error = {out_dir}/{system}.err\n"
        
        "request_memory = 200 GB\n"
        "Queue\n"
    )

    SUBMIT = "condor_submit {job_file}"


def execute(command_string):
    subprocess.Popen(command_string,
                     shell=True,
    )
    
    
def chmod(file_path):
    command = f"chmod 777 {str(file_path)}"
    execute(command)
    
 
def make_dir(dir_path: Path):
    if not dir_path.exists():
        os.mkdir(str(dir_path))
    
    
def write_script(script, file):
    with open(file, "w") as f:
        f.write(script)
    
    
def submit(job_file):
    command = Constants.SUBMIT.format(str(job_file))
    execute(command)
    
    
def main(data_dirs, out_dirs, debug=True):
    
    if debug: print(f"Data dirs: {data_dirs}")
    if debug: print(f"Out dirs: {out_dirs}")
    
    data_dirs_path = Path(data_dirs)
    out_dirs_path = Path(out_dirs)
    
    for data_dir in data_dirs_path.glob("*"):
        if debug: print("#################################################")
        if debug: print(data_dir)
        
        system_name = data_dir.name
        if debug: print(f"System name: {system_name}")
        
        out_dir = out_dirs_path / system_name
        make_dir(out_dir)
        if debug: print(f"Out_dir: {out_dir}")
        
        job_script = Constants.JOB.format(out_dir=str(out_dir), system=system_name)
        if debug: print(f"Job script: {job_script}")
        
        script = Constants.SCRIPT.format(data_dir=str(data_dir), system=system_name)
        if debug: print(f"Script: {script}")
        
        script_file = out_dirs_path / f"{system_name}.sh"
        if debug: print(f"Script file: {script_file}")
        job_file = out_dirs_path / f"{system_name}.job"
        if debug: print(f"job file: {job_file}")
        
        write_script(job_script, job_file)
        write_script(script, script_file)
        
        chmod(script_file)
        
        submit(job_file)        


if __name__ == "__main__":
    fire.Fire(main)
