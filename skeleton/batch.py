from pathlib import Path
import subprocess

import fire

from xlib import database_sql

from constants import Constants


def write(string: str, file: Path):
    with open(file, "w") as f:
        f.write(string)


def submit(command):
    p = subprocess.Popen(
        command,
        shell=True,
    )

    p.communicate()


def dispatch(autobuild: database_sql.Autobuild, out_dir: Path):
    structure = autobuild.path
    event_map = autobuild.event.event_map
    dtag = autobuild.dataset.dtag
    event_idx = autobuild.event.event_idx
    autobuild_id = autobuild.id
    build_id = f"{dtag}_{event_idx}_{autobuild_id}"

    score_script = Constants.SCORE_SCRIPT.format(
        structure=structure,
        event_map=event_map,
        out_dir=str(out_dir),
    )
    score_script_file = out_dir / Constants.SCORE_SCRIPT_FILE.format(build_id=build_id)
    write(score_script, score_script_file)

    executable_file = str(score_script_file)
    log_file = Constants.LOG_FILE.format(build_id=build_id)
    output_file = Constants.OUTPUT_FILE.format(build_id=build_id)
    error_file = Constants.ERROR_FILE.format(build_id=build_id)
    request_memory = Constants.REQUEST_MEMORY
    job_script = Constants.JOB_SCRIPT.format(
        executable_file=executable_file,
        log_file=log_file,
        output_file=output_file,
        error_file=error_file,
        request_memory=request_memory,
    )
    job_script_file = out_dir / Constants.JOB_SCRIPT_FILE
    write(job_script, job_script_file)

    submit_command = Constants.SUBMIT_COMMAND.format()

    submit(submit_command)


def main(database_file: str, out_dir: str):
    database_file = Path(database_file)
    out_dir = Path(out_dir)

    database = database_sql.Database(database_file)

    autobuild_query = database.session.query(database_sql.Autobuild)
    autobuild_list = autobuild_query.all()

    map(
        lambda autobuild: dispatch(autobuild, out_dir),
        autobuild_list,
    )


if __name__ == "__main__":
    fire.Fire(main)
