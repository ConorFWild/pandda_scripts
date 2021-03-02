from __future__ import annotations

from typing import *
from dataclasses import dataclass

from pathlib import Path
import subprocess
import os

import fire

import numpy as np

import gemmi

from constants import Constants

from xlib import database_sql


@dataclass()
class ScoreData:
    pdb: Path
    event_map: Path
    reference: database_sql.ReferenceModel
    event: database_sql.Event


def write(string: str, file: Path):
    with open(file, "w") as f:
        f.write(string)


def submit(command):
    p = subprocess.Popen(
        command,
        shell=True,
    )

    p.communicate()


def chmod(path: Path):
    p = subprocess.Popen(
        f"chmod 777 {str(path)}",
        shell=True,
    )

    p.communicate()


def try_make_dir(path: Path):
    if not path.exists():
        os.mkdir(str(path))


def dispatch(score_data: ScoreData, out_dir: Path):
    structure = score_data.pdb
    event_map = score_data.event_map
    dtag = score_data.reference.dataset.dtag
    event_idx = score_data.event.event_idx
    build_id = f"{dtag}_{event_idx}"

    build_dir = out_dir / build_id

    try_make_dir(build_dir)

    score_script = Constants.SCORE_SCRIPT.format(
        structure=structure,
        event_map=event_map,
        out_dir=str(build_dir),
    )
    score_script_file = build_dir / Constants.SCORE_SCRIPT_FILE.format(build_id=build_id)
    write(score_script, score_script_file)

    chmod(score_script_file)

    executable_file = str(score_script_file)
    log_file = build_dir / Constants.LOG_FILE.format(build_id=build_id)
    output_file = build_dir / Constants.OUTPUT_FILE.format(build_id=build_id)
    error_file = build_dir / Constants.ERROR_FILE.format(build_id=build_id)
    request_memory = Constants.REQUEST_MEMORY
    job_script = Constants.JOB_SCRIPT.format(
        executable_file=executable_file,
        log_file=log_file,
        output_file=output_file,
        error_file=error_file,
        request_memory=request_memory,
    )
    job_script_file = build_dir / Constants.JOB_SCRIPT_FILE.format(build_id=build_id)
    write(job_script, job_script_file)

    submit_command = Constants.SUBMIT_COMMAND.format(job_script_file=str(job_script_file))

    submit(submit_command)


def distance(point_1, point_2):
    dist = np.linalg.norm(np.array(point_1) - np.array(point_2))
    return dist


def get_residue_centroid(residue):
    x = []
    y = []
    z = []

    for atom in residue:
        pos = atom.pos
        x.append(pos.x)
        y.append(pos.y)
        z.append(pos.z)

    centroid = (np.mean(x), np.mean(y), np.mean(z))

    return centroid


def get_ligand_centroid(structure):
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.name == "LIG":
                    centroid = get_residue_centroid(residue)
                    return centroid


def get_scoring_data(database):
    reference_query = database.session.query(database_sql.ReferenceModel)
    event_query = database.session.query(database_sql.Event)

    # Get objects
    reference_list = reference_query.all()
    event_list = event_query.all()

    # Filter to those references with events
    reference_with_events_list = []
    for reference in reference_list:
        reference_events = [event for event in event_list if event.dataset_id == reference.dataset_id]

        if len(reference_events) > 0:
            if Constants.DEBUG:
                print(f"Found {len(reference_events)} events for reference: {reference}")

            reference_with_events_list.append(reference)

        else:
            if Constants.DEBUG:
                print(f"Did not find any events for reference: {reference}")

    if Constants.DEBUG:
        print(
            f"After filtering for references with events, there are {len(reference_with_events_list)} references left to score")

    # Filter to those references with nearby events
    reference_with_nearby_events_list = []
    nearby_event_list = []
    for reference in reference_with_events_list:
        reference_events = [event for event in event_list if event.dataset_id == reference.dataset_id]

        structure = gemmi.read_structure(reference.path)
        ligand_centroid = get_ligand_centroid(structure)
        if Constants.DEBUG:
            print(f"Found reference centroid at : {ligand_centroid}")

        for event in reference_events:

            event_centroid = (event.x, event.y, event.z)
            if Constants.DEBUG:
                print(f"Found event centroid at : {event_centroid}")

            if distance(event_centroid, ligand_centroid) < 3.0:
                reference_with_nearby_events_list.append(reference)
                nearby_event_list.append(event)
                break

    if Constants.DEBUG:
        print(f"After filtering for neraby events, there are {len(nearby_event_list)} references left to score")

    # Construct score datas
    score_data_list = []
    for reference, event in zip(reference_with_nearby_events_list, nearby_event_list):
        pdb = reference.path
        event_map = event.event_map
        score_data = ScoreData(
            pdb,
            event_map,
            reference,
            event,
        )
        score_data_list.append(score_data)

    return score_data_list


def main(database_file: str, out_dir: str):
    database_file = Path(database_file)
    out_dir = Path(out_dir)

    try_make_dir(out_dir)

    database = database_sql.Database(database_file)

    score_data_list: List[ScoreData] = get_scoring_data(database)
    if Constants.DEBUG: print(f"Found {len(score_data_list)} reference event maps to score")

    for score_data in score_data_list:
        if Constants.DEBUG: print(f"Scoring: {score_data.pdb} {score_data.event_map}")
        dispatch(score_data, out_dir)


if __name__ == "__main__":
    fire.Fire(main)
