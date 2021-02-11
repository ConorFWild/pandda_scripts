import os
import subprocess
from pathlib import Path

import fire

from constants import Constants
from xlib.xlib import database_sql


def execute(command: str):
    p = subprocess.Popen(command,
                         shell=True,
                         )

    p.communicate()


def try_make_dir(path: Path):
    if not path.exists():
        os.mkdir(str(path))


def dispatch(event: database_sql.Event, out_dir: Path):
    # Create the Event dir
    event_dir = out_dir / f"{event.dataset.dtag}_{event.event_idx}"
    try_make_dir(event_dir)

    # Write the script that will call python to autobuild the event
    executable_script = Constants.EXECUTABLE.format(event_dir=str(event_dir),
                                                    x=event.x,
                                                    y=event.y,
                                                    z=event.z,
                                                    )
    executable_script_file = event_dir / Constants.EXECUTABLE_SCRIPT_FILE.format(dtag=event.dataset.dtag,
                                                                                 event_idx=event.event_idx)
    with open(executable_script_file, "w") as f:
        f.write(executable_script)

    # Generate a job script file for a condor cluster
    job_script = Constants.JOB.format(executable_script_file=str(executable_script_file))
    job_script_file = event_dir / Constants.JOB_SCRIPT_FILE.format(dtag=event.dataset.dtag,
                                                                   event_idx=event.event_idx)
    with open(job_script_file, "w") as f:
        f.write(job_script)

    # Generate a shell command to submit the job to run the python script
    command = Constants.COMMAND.format(job_script_file=job_script_file)

    # Submit the job
    execute(command)


def main(database_file: str, output_dir: str):
    # Format arguments
    database_file_path = Path(database_file)
    output_dir_path = Path(output_dir)

    # Load database
    database: database_sql.Database = database_sql.Database(database_file_path)

    # Select which datasets to build
    event_query = database.session.query(database_sql.Event)
    event_list = event_query.all()

    # Dispatch to run
    map(
        lambda event: dispatch(event, output_dir_path),
        event_list,
    )


if __name__ == "__main__":
    fire.Fire(main)
