from __future__ import annotations
from re import match

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


# ########
# debug
# #########
def summarise_grid(grid: gemmi.FloatGrid):
    grid_array = np.array(grid, copy=False)
    print((
        f"Grid size: {grid.nu} {grid.nv} {grid.nw} \n"
        f"Grid spacegroup: {grid.spacegroup} \n"
        f"Grid unit cell: {grid.unit_cell} \n"
        f"Grid max: {np.max(grid_array)} \n"
        f"Grid min: {np.min(grid_array)} \n"
        f"Grid mean: {np.mean(grid_array)} \n"
    ))


def summarise_mtz(mtz: gemmi.Mtz):
    mtz_array = np.array(mtz, copy=False)
    print(
        (
            f"Mtz shape: {mtz_array.shape} \n"
            f"Mtz spacegroup: {mtz.spacegroup} \n"
        )
    )


def summarise_structure(structure: gemmi.Structure):
    num_models: int = 0
    num_chains: int = 0
    num_residues: int = 0
    num_atoms: int = 0

    for model in structure:
        num_models += 1
        for chain in model:
            num_chains += 1
            for residue in chain:
                num_residues += 1
                for atom in residue:
                    num_atoms += 1

    print(
        (
            f"Num models: {num_models}"
            f"Num chains: {num_chains}"
            f"Num residues: {num_residues}"
            f"Num atoms: {num_atoms}"
        )
    )


def summarise_event(event: Event):
    print(
        (
            f"Event system: {event.system}\n"
            f"Event dtag: {event.dtag}\n"
            f"Event xyz: {event.x} {event.y} {event.z}\n"
        )
    )


# ########
# data types
# #########

class Constants:
    DEBUG: int = 0
    JOB_SCRIPT = (
        "#################### \n"
        "# \n"
        "# Example 1                                   \n"
        "# Simple HTCondor submit description file \n"
        "#                          \n"
        "####################    \n"

        "Executable   = {executable_file} \n"
        "Log          = {log_file} \n"
        "Output = {output_file} \n"
        "Error = {error_file} \n"

        "request_memory = {request_memory} GB \n"

        "Queue"
    )

    JOB_SCRIPT_FILE = "{system_name}.job"

    LOG_FILE = "{system_name}.log"
    OUTPUT_FILE = "{system_name}.out"
    ERROR_FILE = "{system_name}.err"

    CHANGE_PERMISSION_COMMAND = "chmod 777 {path}"

    SUBMIT_COMMAND = "condor_submit {job_script_file}"

    STRIPPED_RECEPTOR_FILE = "stripped_receptor.pdb"
    LIGAND_FILE = "autobuilding_ligand.cif"
    EVENT_MTZ_FILE = "event.mtz"
    GRAFTED_MTZ_FILE = "grafted.mtz"
    RHOFIT_DIR = "rhofit"
    RHOFIT_EVENT_DIR = "rhofit_{}"
    RHOFIT_NORMAL_DIR = "rhofit_normal"
    RHOFIT_RESULTS_FILE = "results.txt"
    RHOFIT_RESULT_JSON_FILE = "result.json"
    RHOFIT_NORMAL_RESULT_JSON_FILE = "result_normal.json"
    RHOFIT_BEST_MODEL_FILE = "best.pdb"
    RSCC_TABLE_FILE = "rscc_table.csv"

    BUILD_DIR_PATTERN = "{pandda_name}_{dtag}_{event_idx}"

    CUMULATIVE_HITS_PLOT_FILE = "cumulative_hits.png"

    PANDDA_ANALYSES_DIR = "analyses"
    PANDDA_ANALYSE_EVENTS_FILE = "pandda_analyse_events.csv"
    PANDDA_ANALYSE_SITES_FILE = "pandda_analyse_sites.csv"
    PANDDA_PROCESSED_DATASETS_DIR = "processed_datasets"
    PANDDA_MODELLED_STRUCTURES_DIR = "modelled_structures"
    PANDDA_LIGAND_FILES_DIR = "ligand_files"
    PANDDA_PDB_FILE = "{}-pandda-input.pdb"
    PANDDA_MTZ_FILE = "{}-pandda-input.mtz"

    PANDDA_LIGAND_CIF_FILE = "ligand.cif"
    PANDDA_LIGAND_PDB_FILE = "ligand.pdb"
    PANDDA_LIGAND_SMILES_FILE = "ligand.smiles"

    PANDDA_INSPECT_EVENTS_PATH = "pandda_inspect_events.csv"
    PANDDA_EVENT_MAP_FILE = "{}-event_{}_1-BDC_{}_map.native.ccp4"
    PANDDA_EVENT_MODEL = "{}-pandda-model.pdb"

    PANDDA_PROTEIN_MASK_FILE = "protein_mask.ccp4"
    PANDDA_SYMMETRY_MASK_FILE = "symmetry_mask.ccp4"
    PANDDA_MEAN_MAP_FILE = "mean_{number}_{res}.ccp4"
    PANDDA_SIGMA_S_M_FILE = "sigma_s_m_{number}_{res}.ccp4"

    PANDDA_Z_MAP_FILE = "{dtag}-z_map.native.ccp4"
    PANDDA_EVENT_MAP_FILE = "{dtag}-event_{event_idx}_1-BDC_{bdc}_map.native.ccp4"

    PANDDA_LOG_FILE = "pandda_log.json"

    ELBOW_COMMAND = "cd {event_dir.event_dir}; phenix. {smiles_file.smiles_file} --output=\"{prefix}\""
    FINAL_LIGAND_CIF_FILE = "final_ligand_cif.cif"

    MASKED_PDB_FILE = "masked.pdb"

    CUT_OUT_EVENT_MAP_FILE = "cut_out_event_map.ccp4"

    RHOFIT_LOG_FILE = "rhofit.log"
    RHOFIT_OUTPUT_FILE = "rhofit.out"
    RHOFIT_ERROR_FILE = "rhofit.err"
    RHOFIT_JOB_SCRIPT_FILE = "rhofit.job"

    RHOFIT_SCRIPT = (
        "#!/bin/bash \n"
        "conda activate pandda \n"
        "source ~/.bashrc \n"
        "source /data/share-2/conor/xtal_software/ccp4-7.1/bin/ccp4.setup-sh \n"
        "source /data/share-2/conor/xtal_software/phenix/phenix-1.18.2-3874/phenix_env.sh \n"
        "source /data/share-2/conor/xtal_software/buster-2.10/setup.sh \n"
        "rhofit -m {mtz} -p {pdb} -l {ligand} -d {out_dir_path} -allclusters -use_2fofc -thorough "
    )
    RHOFIT_SCRIPT_CLEMENS = (
        "#!/bin/bash \n"
        "conda activate pandda \n"
        "source ~/.bashrc \n"
        "source /data/share-2/conor/xtal_software/ccp4-7.1/bin/ccp4.setup-sh \n"
        "source /data/share-2/conor/xtal_software/phenix/phenix-1.18.2-3874/phenix_env.sh \n"
        "source /data/share-2/conor/xtal_software/buster-2.10/setup.sh \n"
        "{pandda_rhofit} -map {event_map} -mtz {mtz} -pdb {pdb} -cif {cif} -out {out_dir_path}"
    )
    PANDDA_RHOFIT_SCRIPT_FILE = "/data/share-2/conor/pandda/pandda_scripts/pandda_rhofit.sh"
    RHOFIT_SCRIPT_FILE = "run_rhofit.sh"

    PHASE_GRAFTED_MTZ_FILE = "phase_grafted_mtz.mtz"

    RHOFIT_HIT_REGEX = "(Hit_[^\s]+)[\s]+[^\s]+[\s]+[^\s]+[\s]+([^\s]+)"
    RHOFIT_CLUSTER_BUILD_REGEX = "Hit_([^_]+)_[^_]+_([^_]+).pdb"

    RHOFIT_EVENT_BUILD_RESULT_FILE = "result.json"

    CLUSTER = "CONDOR"
    CONDOR_JOB_ID_REGEX = r"[0-9]+\."
    CONDOR_STATUS_COMMAND = r"condor_q"


@dataclasses.dataclass()
class Args:
    pandda_dirs_dir: Path
    autobuild_dirs_dir: Path
    reference_structure_dir: Path
    debug: int

    @staticmethod
    def from_args(args: Any):
        pandda_dirs_dir = Path(args.pandda_dirs_dir)
        autobuild_dirs_dir = Path(args.autobuild_dirs_dir)
        reference_structure_dir = Path(args.reference_structure_dir)
        debug: int = int(args.debug)

        return Args(
            pandda_dirs_dir,
            autobuild_dirs_dir,
            reference_structure_dir,
            debug,
        )

    @staticmethod
    def from_cmd():
        parser = argparse.ArgumentParser()
        parser.add_argument("--pandda_dirs_dir",
                            )
        parser.add_argument("--autobuild_dirs_dir",
                            )
        parser.add_argument("--reference_structure_dir",
                            )
        parser.add_argument("--debug",
                            default=2,
                            )
        args = parser.parse_args()

        return Args.from_args(args)


@dataclasses.dataclass()
class HKL:
    h: int
    k: int
    l: int

    @staticmethod
    def from_hkl(h, k, l):
        return HKL(h, k, l)

    def __hash__(self):
        return hash((self.h, self.k, self.l))

    def to_list(self) -> List[int]:
        return [self.h, self.k, self.l]

    @staticmethod
    def from_list(hkl_list: List[int]):
        return HKL(hkl_list[0],
                   hkl_list[1],
                   hkl_list[2],
                   )

    def is_000(self):
        if self.h == 0:
            if self.k == 0:
                if self.l == 0:
                    return True

        return False


@dataclasses.dataclass()
class Reflection:
    hkl: HKL
    data: np.array

    @staticmethod
    def from_row(row: np.array):
        h = int(row[0])
        k = int(row[1])
        l = int(row[2])
        hkl = HKL(h, k, l)
        data = row[3:]
        return Reflection(hkl,
                          data)


@dataclasses.dataclass()
class ReflectionsDict:
    reflections_dict: Dict[HKL, Reflection]

    @staticmethod
    def from_array(array):
        reflections = {}
        for row in array:
            reflection = Reflection.from_row(row)
            reflections[reflection.hkl] = reflection

        return ReflectionsDict(reflections)

    def __getitem__(self, item):
        return self.reflections_dict[item]

    def __iter__(self):
        for hkl in self.reflections_dict:
            yield hkl

    def to_array(self):
        rows = []
        for hkl in self.reflections_dict:
            hkl_array = np.array(hkl.to_list())
            data = self.reflections_dict[hkl].data

            row = np.hstack([hkl_array, data])

            rows.append(row)

        array = np.vstack(rows)

        return array


@dataclasses.dataclass()
class System:
    system: str

    def __hash__(self) -> int:
        return hash(self.system)


@dataclasses.dataclass()
class Dtag:
    dtag: str

    def __hash__(self) -> int:
        return hash(self.dtag)


@dataclasses.dataclass()
class EventIDX:
    event_idx: str

    def __hash__(self) -> int:
        return hash(self.event_idx)


@dataclasses.dataclass()
class EventID:
    system: System
    dtag: Dtag
    event_idx: EventIDX

    def __hash__(self) -> int:
        return hash(
            (self.system,
             self.dtag,
             self.event_idx,
             )
        )


@dataclasses.dataclass()
class Event:
    event_input_dir: Path
    event_output_dir: Path
    system: System
    dtag: Dtag
    event_idx: EventIDX
    bdc: float
    x: float
    y: float
    z: float


@dataclasses.dataclass()
class BuildNumberID:
    build_number_id: int

    def __hash__(self):
        return hash(self.build_number_id)


@dataclasses.dataclass()
class BuildClusterID:
    build_cluster_id: int

    def __hash__(self):
        return hash(self.build_cluster_id)


@dataclasses.dataclass()
class Build:
    build_file: Path
    build_rscc: float


@dataclasses.dataclass()
class ClusterBuildResults:
    cluster_build_results: Dict[BuildNumberID, Build]

    def __iter__(self):
        for build_number_id in self.cluster_build_results:
            yield build_number_id

    def __getitem__(self, item):
        return self.cluster_build_results[item]


@dataclasses.dataclass()
class EventBuildResults:
    success: bool
    build_results: Dict[BuildClusterID, ClusterBuildResults]

    def __iter__(self):
        for build_cluster_id in self.build_results:
            yield build_cluster_id

    def __getitem__(self, item):
        return self.build_results[item]

    @classmethod
    def from_event(cls, event: Event):

        rhofit_dir: Path = event.event_output_dir / Constants.RHOFIT_DIR

        rhofit_results_file: Path = rhofit_dir / Constants.RHOFIT_RESULTS_FILE

        if not rhofit_results_file.exists():
            return EventBuildResults(False, {})

        with open(str(rhofit_results_file), "r") as f:
            results_string = f.read()

        build_matches = re.findall(Constants.RHOFIT_HIT_REGEX,
                                   results_string,
                                   )

        cluster_builds = {}
        for build_match in build_matches:
            build_file = build_match[0]
            rscc = build_match[1]

            cluster_build_matches = re.findall(Constants.RHOFIT_CLUSTER_BUILD_REGEX,
                                               build_file,
                                               )
            cluster = BuildClusterID(int(cluster_build_matches[0][0]))
            build_number = BuildNumberID(int(cluster_build_matches[0][1]))

            if cluster not in cluster_builds:
                cluster_builds[cluster] = {}

            cluster_builds[cluster][build_number] = Build(build_file=rhofit_dir / build_file,
                                                          build_rscc=rscc,
                                                          )

        return EventBuildResults.from_dict(cluster_builds)

    @classmethod
    def from_dict(cls, cluster_builds_dict: Dict[BuildClusterID, Dict[BuildNumberID, Build]]):
        event_clusters = {}
        for cluster_id in cluster_builds_dict:
            cluster_builds = cluster_builds_dict[cluster_id]
            event_clusters[cluster_id] = ClusterBuildResults(cluster_builds)

        return EventBuildResults(True, event_clusters)

    def to_json_file(self, event: Event):
        builds = {}
        clusters = self.build_results
        for cluster_id in clusters:
            if cluster_id.build_cluster_id not in builds:
                builds[cluster_id.build_cluster_id] = {}

            build_results = clusters[cluster_id]
            for build_number in build_results:
                builds[cluster_id.build_cluster_id][
                    build_number.build_number_id] = {
                    "build_file": str(build_results[build_number].build_file),
                    "build_rscc": float(build_results[build_number].build_rscc),
                }

        result_json_file: Path = event.event_output_dir / Constants.RHOFIT_RESULT_JSON_FILE

        with open(result_json_file, "w") as f:
            json.dump(builds, f)

        return result_json_file


# ########
# Get Files functions
# #########

def get_event_input_dir(system: System, dtag: Dtag, event_idx: EventIDX, pandda_dirs_dir: Path) -> Path:
    return pandda_dirs_dir / system.system / Constants.PANDDA_PROCESSED_DATASETS_DIR / dtag.dtag


def get_pdb_file(event: Event) -> Path:
    event_dir: Path = event.event_input_dir
    pdb_file: Path = event_dir / Constants.PANDDA_PDB_FILE.format(event.dtag.dtag)
    return pdb_file


def get_ligand_cif_file(event: Event) -> Union[Path, None]:
    event_dir: Path = event.event_input_dir
    ligand_cif_file = event_dir / Constants.PANDDA_LIGAND_CIF_FILE
    if ligand_cif_file.exists():
        return ligand_cif_file
    else:
        return None


def get_ligand_smiles_file(event: Event) -> Union[Path, None]:
    event_dir: Path = event.event_input_dir
    compound_dir: Path = event_dir / Constants.PANDDA_LIGAND_FILES_DIR
    ligand_smiles_file = compound_dir / Constants.PANDDA_LIGAND_SMILES_FILE
    if ligand_smiles_file.exists():
        return ligand_smiles_file
    else:
        return None


def get_ligand_pdb_file(event: Event) -> Union[Path, None]:
    event_dir: Path = event.event_input_dir
    compound_dir: Path = event_dir / Constants.PANDDA_LIGAND_FILES_DIR
    ligand_pdb_file = compound_dir / Constants.PANDDA_LIGAND_PDB_FILE
    if ligand_pdb_file.exists():
        return ligand_pdb_file
    else:
        return None


def get_mtz_file(event: Event) -> Path:
    event_dir: Path = event.event_input_dir
    mtz_file: Path = event_dir / Constants.PANDDA_MTZ_FILE.format(event.dtag.dtag)
    return mtz_file


def get_event_map_file(event: Event) -> Path:
    event_dir: Path = event.event_input_dir
    event_map_file: Path = event_dir / Constants.PANDDA_EVENT_MAP_FILE.format(
        dtag=event.dtag.dtag,
        event_idx=event.event_idx.event_idx,
        bdc=event.bdc,
    )
    return event_map_file


# ########
# Make inpuit files
# #########

def get_event_output_dir(system: System, dtag: Dtag, event_idx: EventIDX, autobuild_dirs_dir: Path) -> Path:
    return autobuild_dirs_dir / system.system / dtag.dtag / event_idx.event_idx


def get_pdb(pdb_file: Path) -> gemmi.Structure:
    structure: gemmi.Structure = gemmi.read_structure(str(pdb_file))
    return structure


def get_event_map(event_map_file: Path) -> gemmi.FloatGrid:
    m = gemmi.read_ccp4_map(str(event_map_file))
    m.setup()

    grid_array = np.array(m.grid, copy=True)

    new_grid = gemmi.FloatGrid(*grid_array.shape)
    new_grid.spacegroup = m.grid.spacegroup  # gemmi.find_spacegroup_by_name("P 1")
    new_grid.set_unit_cell(m.grid.unit_cell)

    new_grid_array = np.array(new_grid, copy=False)
    new_grid_array[:, :, :] = grid_array[:, :, :]

    return new_grid


def get_mtz(mtz_file: Path) -> gemmi.Mtz:
    mtz: gemmi.Mtz = gemmi.read_mtz_file(str(mtz_file))
    return mtz


def get_masked_pdb(pdb: gemmi.Structure, event: Event, radius: float = 7.0) -> gemmi.Structure:
    event_centoid = gemmi.Position(
        event.x,
        event.y,
        event.z,
    )

    new_structure = gemmi.Structure()

    for model_i, model in enumerate(pdb):
        new_model = gemmi.Model(model.name)
        new_structure.add_model(new_model, pos=-1)

        for chain_i, chain in enumerate(model):
            new_chain = gemmi.Chain(chain.name)
            new_structure[model_i].add_chain(new_chain, pos=-1)

            for residue_i, residue in enumerate(chain):
                new_residue = gemmi.Residue()
                new_residue.name = residue.name
                new_residue.seqid = residue.seqid
                new_residue.subchain = residue.subchain
                new_residue.label_seq = residue.label_seq
                new_residue.het_flag = residue.het_flag
                new_structure[model_i][chain_i].add_residue(new_residue, pos=-1)

                for atom_i, atom in enumerate(residue):
                    pos = atom.pos
                    if pos.dist(event_centoid) > radius:
                        new_structure[model_i][chain_i][residue_i].add_atom(atom, pos=-1)

    for model_i, model in enumerate(pdb):
        pdb.add_model(new_structure[model_i], pos=-1)
        del pdb[0]

    return pdb


def get_masked_pdb_file(masked_pdb: gemmi.Structure, event: Event) -> Path:
    masked_pdb_file: Path = event.event_output_dir / Constants.MASKED_PDB_FILE
    masked_pdb.write_minimal_pdb(str(masked_pdb_file))

    return masked_pdb_file


def get_cut_out_event_map(event: Event, event_map: gemmi.FloatGrid, radius: float = 10.0) -> gemmi.FloatGrid:
    event_centroid = gemmi.Position(event.x, event.y, event.z)

    xmap_array = np.array(event_map, copy=True)

    mask_grid = gemmi.Int8Grid(*xmap_array.shape)
    mask_grid.spacegroup = event_map.spacegroup
    mask_grid.set_unit_cell(event_map.unit_cell)

    mask_grid.set_points_around(event_centroid,
                                radius=radius,
                                value=1,
                                )
    mask_grid.symmetrize_max()

    mask_array = np.array(mask_grid, copy=False, dtype=np.int8)

    new_grid = gemmi.FloatGrid(*xmap_array.shape)
    new_grid.spacegroup = event_map.spacegroup  # gemmi.find_spacegroup_by_name("P 1")
    new_grid.set_unit_cell(event_map.unit_cell)

    new_grid_array = np.array(new_grid, copy=False)

    new_grid_array[np.nonzero(mask_array)] = xmap_array[np.nonzero(mask_array)]
    new_grid.symmetrize_max()

    return new_grid


def get_cut_out_event_mtz(cut_out_event_map: gemmi.FloatGrid, mtz: gemmi.Mtz) -> gemmi.Mtz:
    # cut_out_event_map.spacegroup = mtz.spacegroup
    cut_out_event_map.spacegroup = mtz.spacegroup

    cut_out_event_map.symmetrize_max()

    sf = gemmi.transform_map_to_f_phi(cut_out_event_map, half_l=False)
    data = sf.prepare_asu_data(dmin=mtz.resolution_high(), with_000=True)

    cut_out_event_mtz = gemmi.Mtz(with_base=True)
    cut_out_event_mtz.spacegroup = sf.spacegroup
    cut_out_event_mtz.cell = sf.unit_cell
    cut_out_event_mtz.add_dataset('unknown')
    cut_out_event_mtz.add_column('FWT', 'F')
    cut_out_event_mtz.add_column('PHWT', 'P')
    cut_out_event_mtz.set_data(data)

    return cut_out_event_mtz


def phase_graft(mtz: gemmi.Mtz, cut_out_event_mtz: gemmi.Mtz) -> gemmi.Mtz:
    initial_mtz = mtz
    event_mtz = cut_out_event_mtz

    initial_mtz_data = np.array(initial_mtz, copy=False)
    event_mtz_data = np.array(event_mtz, copy=False)

    initial_reflections: ReflectionsDict = ReflectionsDict.from_array(initial_mtz_data)
    event_reflections: ReflectionsDict = ReflectionsDict.from_array(event_mtz_data)

    initial_asu = gemmi.ReciprocalAsu(initial_mtz.spacegroup)
    operations = initial_mtz.spacegroup.operations()

    initial_mtz_fwt_index = initial_mtz.column_labels().index("FWT")
    event_mtz_fwt_index = event_mtz.column_labels().index("FWT")

    initial_mtz_phwt_index = initial_mtz.column_labels().index("PHWT")
    event_mtz_phwt_index = event_mtz.column_labels().index("PHWT")

    fom_index = initial_mtz.column_labels().index("FOM")
    initial_mtz_fo_index = initial_mtz.column_labels().index("F")
    initial_mtz_fc_index = initial_mtz.column_labels().index("FC")
    initial_mtz_phc_index = initial_mtz.column_labels().index("PHIC")
    initial_mtz_r_index = initial_mtz.column_labels().index("FreeR_flag")
    initial_mtz_ls_fc_all_index = initial_mtz.column_labels().index("FC_ALL_LS")
    initial_mtz_ls_phc_all_index = initial_mtz.column_labels().index("PHIC_ALL_LS")
    initial_mtz_fc_all_index = initial_mtz.column_labels().index("FC_ALL")
    initial_mtz_phc_all_index = initial_mtz.column_labels().index("PHIC_ALL")
    initial_mtz_delfwt_index = initial_mtz.column_labels().index("DELFWT")
    initial_mtz_phdelwt_index = initial_mtz.column_labels().index("PHDELWT")

    initial_mtz_sigf_index = initial_mtz.column_labels().index("SIGF")

    print("\tBeginning graft...")
    new_reflections = {}
    for hkl in event_reflections:
        event_reflection = event_reflections[hkl]

        asu_hkl_list, asu_hkl_value = initial_asu.to_asu(hkl.to_list(), operations, )

        asu_hkl = HKL.from_list(asu_hkl_list)
        if asu_hkl.is_000():
            data = np.zeros(len(list(initial_reflections.reflections_dict.values())[0].data))
            new_reflection = Reflection(hkl, data)

        # elif asu_hkl not in initial_reflections:
        #     print(f"\tMissing reflection: {asu_hkl}")
        #     continue

        else:
            try:
                initial_reflection: Reflection = initial_reflections[asu_hkl]

                new_reflection = Reflection(hkl, np.copy(initial_reflection.data))
            except Exception as e:
                print(f"\tWARNING: Missing Reflection: {asu_hkl}")
                data = np.zeros(len(list(initial_reflections.reflections_dict.values())[0].data))
                new_reflection = Reflection(hkl, data)

        new_reflection.data[initial_mtz_fwt_index - 3] = event_reflection.data[event_mtz_fwt_index - 3]
        new_reflection.data[initial_mtz_phwt_index - 3] = event_reflection.data[event_mtz_phwt_index - 3]

        new_reflections[hkl] = new_reflection

    new_array = ReflectionsDict(new_reflections).to_array()

    initial_mtz.spacegroup = event_mtz.spacegroup
    initial_mtz.set_data(new_array)

    return initial_mtz


def get_cut_out_event_map_file(event: Event, get_cut_out_event_map: gemmi.FloatGrid, p1: bool = True) -> Path:
    cut_out_event_file: Path = event.event_output_dir / Constants.CUT_OUT_EVENT_MAP_FILE
    ccp4 = gemmi.Ccp4Map()
    ccp4.grid = get_cut_out_event_map
    if p1:
        ccp4.grid.spacegroup = gemmi.find_spacegroup_by_name("P 1")
    else:
        ccp4.grid.symmetrize_max()
    ccp4.update_ccp4_header(2, True)
    ccp4.write_ccp4_map(str(cut_out_event_file))

    return cut_out_event_file


def write_mtz(mtz: gemmi.Mtz, path: Path):
    mtz.write_to_file(str(path))


def get_phase_grafted_mtz_file(event: Event, phase_grafted_mtz: gemmi.Mtz) -> Path:
    phase_grafted_mtz_file: Path = event.event_output_dir / Constants.PHASE_GRAFTED_MTZ_FILE
    phase_grafted_mtz.write_to_file(str(phase_grafted_mtz_file))
    return phase_grafted_mtz_file


def get_elbow_command(input_file: Path) -> str:
    command = Constants.ELBOW_COMMAND
    return command


def submit(command: str) -> Tuple[str, str]:
    p = subprocess.Popen(
        command,
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    stdin, stderr = p.communicate()

    return str(stdin), str(stderr)


def get_final_ligand_cif_file(
        event: Event,
        ligand_cif_file: Union[Path, None],
        ligand_smiles_file: Union[Path, None],
        ligand_pdb_file: Union[Path, None],
) -> Union[Path, None]:
    event_output_dir: Path = event.event_output_dir
    final_ligand_cif_file: Path = event_output_dir / Constants.FINAL_LIGAND_CIF_FILE

    if ligand_pdb_file:
        command: str = get_elbow_command(ligand_pdb_file)
        stdout, stderr = submit(command)
        return final_ligand_cif_file

    elif ligand_smiles_file:
        command: str = get_elbow_command(ligand_smiles_file)
        stdout, stderr = submit(command)
        return final_ligand_cif_file


    elif ligand_cif_file:
        final_ligand_cif_file: Path = cast(Path, ligand_cif_file)
        return Path(ligand_cif_file)

    else:
        # raise Exception("No ligand!")
        return None


def change_permission(path: Path):
    command: str = Constants.CHANGE_PERMISSION_COMMAND.format(path=path)
    submit(command)


# ########
# Get run rhofit
# #########

def get_rhofit_script_clemens(
        event: Event,
        event_map_file: Path,
        mtz_file: Path,
        masked_pdb_file: Path,
        ligand_cif_file: Path,
) -> str:
    event_output_dir: Path = event.event_output_dir
    # rhofit_dir: Path = event_output_dir / Constants.RHOFIT_DIR 
    out_dir: Path = event_output_dir
    rhofit_command: str = Constants.RHOFIT_SCRIPT_CLEMENS.format(
        pandda_rhofit=Constants.PANDDA_RHOFIT_SCRIPT_FILE,
        event_map=event_map_file,
        mtz=mtz_file,
        pdb=masked_pdb_file,
        cif=ligand_cif_file,
        out_dir_path=out_dir,
    )
    return rhofit_command


def get_rhofit_script_file(event: Event, rhofit_script: str) -> Path:
    rhofit_script_file: Path = event.event_output_dir / Constants.RHOFIT_SCRIPT_FILE
    with open(str(rhofit_script_file), "w") as f:
        f.write(rhofit_script)

    return rhofit_script_file


def get_job_script(event: Event, rhofit_script_file: Path, request_memory: int = 15) -> str:
    event_output_dir: Path = event.event_output_dir
    executable_file: Path = rhofit_script_file
    log_file = event_output_dir / Constants.RHOFIT_LOG_FILE
    output_file = event_output_dir / Constants.RHOFIT_OUTPUT_FILE
    error_file = event_output_dir / Constants.RHOFIT_ERROR_FILE

    job_script: str = Constants.JOB_SCRIPT.format(
        executable_file=executable_file,
        log_file=log_file,
        output_file=output_file,
        error_file=error_file,
        request_memory=request_memory,
    )

    return job_script


def get_job_script_file(job_script: str, event: Event) -> Path:
    event_output_dir: Path = event.event_output_dir
    job_script_file: Path = event_output_dir / Constants.RHOFIT_JOB_SCRIPT_FILE
    with open(job_script_file, "w") as f:
        f.write(job_script)

    return job_script_file


def get_submit_command(job_script_file: Path) -> str:
    submit_command: str = Constants.SUBMIT_COMMAND.format(job_script_file=job_script_file)
    return submit_command


def get_condor_jobid(stdout: str) -> str:
    pattern: Pattern = re.compile(Constants.CONDOR_JOB_ID_REGEX)
    job_id_match_list: List[str] = pattern.findall(stdout)
    return job_id_match_list[0]


def get_condor_cluster_status() -> str:
    command: str = Constants.CONDOR_STATUS_COMMAND
    stdout, stderr = submit(command)

    return stdout


def get_job_running(jobid_str: str, cluster_status_str: str) -> bool:
    pattern: Pattern = re.compile(jobid_str)
    matches: List[str] = pattern.findall(cluster_status_str)

    if len(matches) == 0:
        return False

    else:
        return True


def check_job(stdout: str):
    if Constants.CLUSTER == "CONDOR":
        jobid_str: str = get_condor_jobid(stdout)

    if Constants.DEBUG > 0: print(f"Jobid is: {jobid_str}")

    while True:
        if Constants.CLUSTER == "CONDOR":
            cluster_status_str: str = get_condor_cluster_status()

        if Constants.DEBUG > 0: print(f"cluster_status_str is: {cluster_status_str}")

        if get_job_running(jobid_str, cluster_status_str):
            if Constants.DEBUG > 0: print(f"\tNot yet finished!")
            time.sleep(5)
        else:
            break


# ########
# Script functions
# #########

def build_event(event: Event):
    # ########
    # Debug event info
    # #########    
    if Constants.DEBUG > 0: summarise_event(event)

    # ########
    # Check if there is already an event
    # #########    
    build_result_file: Path = event.event_output_dir / Constants.RHOFIT_RESULT_JSON_FILE
    if Constants.DEBUG > 0: print(f"Build result file is {build_result_file}")
    if build_result_file.exists():
        with open(str(build_result_file), "r") as f:
            result_dict = json.load(f)

        return EventBuildResults.from_dict(result_dict)

    # ########
    # Get Files
    # #########

    # Get pdb path
    pdb_file: Path = get_pdb_file(event)
    if Constants.DEBUG > 0: print(f"pdb_file is: {pdb_file}")

    # Get ligand cifs paths
    ligand_cif_file: Union[Path, None] = get_ligand_cif_file(event)
    if Constants.DEBUG > 0: print(f"ligand_cif_file is: {ligand_cif_file}")

    # Get ligand smile paths
    ligand_smiles_file: Union[Path, None] = get_ligand_smiles_file(event)
    if Constants.DEBUG > 0: print(f"ligand_smiles_file is: {ligand_smiles_file}")

    # Get ligand pdb paths
    ligand_pdb_file: Union[Path, None] = get_ligand_pdb_file(event)
    if Constants.DEBUG > 0: print(f"ligand_pdb_file is: {ligand_pdb_file}")

    # Get original mtz paths
    mtz_file: Path = get_mtz_file(event)
    if Constants.DEBUG > 0: print(f"mtz_file is: {mtz_file}")

    # Event map file
    event_map_file: Path = get_event_map_file(event)
    if Constants.DEBUG > 0: print(f"event_map_file is: {event_map_file}")

    # copy over if debug
    if Constants.DEBUG > 0:
        initial_pdb_file: Path = event.event_output_dir / "initial_pdb.pdb"
        initial_mtz_file: Path = event.event_output_dir / "initial_mtz.mtz"
        initial_event_map_file: Path = event.event_output_dir / "event_map.ccp4"

        shutil.copyfile(pdb_file, initial_pdb_file)
        shutil.copyfile(mtz_file, initial_mtz_file)
        shutil.copyfile(event_map_file, initial_event_map_file)

    # ########
    # Make input files
    # #########

    # Get pdb
    pdb: gemmi.Structure = get_pdb(pdb_file)
    if Constants.DEBUG > 0:
        print("# initial pdb")
        summarise_structure(pdb)

    # Get event map 
    event_map: gemmi.FloatGrid = get_event_map(event_map_file)
    if Constants.DEBUG > 0:
        print("# Event map")
        summarise_grid(event_map)

    # Cut out events:
    cut_out_event_map: gemmi.FloatGrid = get_cut_out_event_map(event, event_map)
    if Constants.DEBUG > 0:
        print("# Cut event map")
        summarise_grid(cut_out_event_map)

    # Save cut out event
    cut_out_event_map_file: Path = get_cut_out_event_map_file(event, cut_out_event_map)

    # Get ligand file
    final_ligand_cif_file: Union[Path, None] = get_final_ligand_cif_file(
        event,
        ligand_cif_file,
        ligand_smiles_file,
        ligand_pdb_file,
    )

    if final_ligand_cif_file is None:
        return None
    else:
        cast(Path, final_ligand_cif_file)

    # ########
    # Run Rhofit
    # #########

    # Make rhofit commands
    rhofit_script: str = get_rhofit_script_clemens(
        event,
        cut_out_event_map_file,
        mtz_file,
        pdb_file,
        final_ligand_cif_file,
    )
    if Constants.DEBUG > 0: print(f"rhofit_command is: {rhofit_script}")

    # Rhofit script file
    rhofit_script_file: Path = get_rhofit_script_file(event, rhofit_script)
    if Constants.DEBUG > 0: print(f"rhofit_script_file is: {rhofit_script_file}")

    # Change permissions
    change_permission(rhofit_script_file)

    # Make job script
    job_script: str = get_job_script(event, rhofit_script_file)
    if Constants.DEBUG > 0: print(f"job_script is: {job_script}")

    # Write job script
    job_script_file: Path = get_job_script_file(job_script, event)
    if Constants.DEBUG > 0: print(f"job_script_file is: {job_script_file}")

    # Submit command
    submit_command: str = get_submit_command(job_script_file)
    if Constants.DEBUG > 0: print(f"submit_command is: {submit_command}")

    # Execute job script
    stdout, stderr = submit(submit_command)
    if Constants.DEBUG > 0: print(f"stdin is: {stdout}")
    if Constants.DEBUG > 0: print(f"stderr is: {stderr}")

    # Check if job is finished
    check_job(stdout)

    # Parse Rhofit result
    result: EventBuildResults = EventBuildResults.from_event(event)
    if Constants.DEBUG > 0: print(f"result is: {result}")

    # Save result
    result_json_file: Path = result.to_json_file(event)
    if Constants.DEBUG > 0: print(f"result_json_file is: {result_json_file}")

    # return
    return result_json_file


def get_event_table_dict(path_list: List[Path]) -> Dict[System, pd.DataFrame]:
    event_table_dict = {}

    for path in path_list:
        system: System = System(path.name)

        event_table_file: Path = path / Constants.PANDDA_ANALYSES_DIR / Constants.PANDDA_ANALYSE_EVENTS_FILE

        if Constants.DEBUG > 0: print(event_table_file)

        if event_table_file.exists():
            # try:
            event_table: pd.DataFrame = pd.read_csv(str(event_table_file))
            event_table_dict[system] = event_table
            # except Exception as e:
            #     print(e)
            #     print(f"event_table_file seems empty?")
            #     continue

    return event_table_dict


def get_event_id(system: System, row: pd.Series) -> EventID:
    dtag: Dtag = Dtag(row["dtag"])
    event_idx: EventIDX = EventIDX(row["event_idx"])
    return EventID(
        system=system,
        dtag=dtag,
        event_idx=event_idx,
    )


def get_event(system: System, row: pd.Series, pandda_dirs_dir: Path, autobuild_dirs_dir: Path) -> Event:
    dtag: Dtag = Dtag(row["dtag"])
    event_idx: EventIDX = EventIDX(str(row["event_idx"]))
    bdc = row["1-BDC"]
    x = row["x"]
    y = row["y"]
    z = row["z"]

    event_input_dir: Path = get_event_input_dir(system, dtag, event_idx, pandda_dirs_dir)
    event_output_dir: Path = get_event_output_dir(system, dtag, event_idx, autobuild_dirs_dir)

    return Event(
        event_input_dir,
        event_output_dir,
        system=system,
        dtag=dtag,
        event_idx=event_idx,
        bdc=bdc,
        x=x,
        y=y,
        z=z,
    )


def get_event_dict(event_table_dict: Dict[System, pd.DataFrame], pandda_dirs_dir: Path, autobuild_dirs_dir: Path) -> \
Dict[EventID, Event]:
    event_dict: Dict[EventID, Event] = {}

    for system, event_table in event_table_dict.items():
        for index, row in event_table.iterrows():
            event_id: EventID = get_event_id(system, row)
            event: Event = get_event(system, row, pandda_dirs_dir, autobuild_dirs_dir)

            event_dict[event_id] = event

    return event_dict


def make_system_dir(autobuild_dirs_dir: Path, event: Event):
    system_dir = autobuild_dirs_dir / event.system.system

    if system_dir.exists():
        return

    else:
        mkdir(str(system_dir))


def make_dtag_dir(autobuild_dirs_dir: Path, event: Event):
    dtag_dir = autobuild_dirs_dir / event.system.system / event.dtag.dtag

    if dtag_dir.exists():
        return

    else:
        mkdir(str(dtag_dir))


def make_event_idx_dir(autobuild_dirs_dir: Path, event: Event):
    event_idx_dir = autobuild_dirs_dir / event.system.system / event.dtag.dtag / event.event_idx.event_idx

    if event_idx_dir.exists():
        return

    else:
        mkdir(str(event_idx_dir))


def make_rhofit_dir(autobuild_dirs_dir: Path, event: Event):
    event_rhofit_dir: Path = autobuild_dirs_dir / event.system.system / event.dtag.dtag / event.event_idx.event_idx / Constants.RHOFIT_DIR

    if event_rhofit_dir.exists():
        return

    else:
        mkdir(str(event_rhofit_dir))


def make_event_output_dir(event: Event, autobuild_dirs_dir: Path):
    make_system_dir(autobuild_dirs_dir, event)
    make_dtag_dir(autobuild_dirs_dir, event)
    make_event_idx_dir(autobuild_dirs_dir, event)
    # make_rhofit_dir(autobuild_dirs_dir, event)


def make_autobuild_output_dir(event_dict: Dict[EventID, Event], autobuild_dirs_dir: Path):
    for event_id, event in event_dict.items():
        make_event_output_dir(event, autobuild_dirs_dir)


def map_dict(func: Callable, dictionary: Dict[EventID, Event]):
    # if Constants.DEBUG > 0:
    #     for value in dictionary.values():
    #         func(value)
    # else:
    values = list(dictionary.values())
    joblib.Parallel(
        n_jobs=1,
        verbose=50,
    )(
        joblib.delayed(
            func
        )(
            value
        )
        for value
        in values
    )


def get_event_with_reference_dict(
        event_dict,
        reference_structure_dict: ReferenceStructureDict,
) -> Dict[EventID, Event]:
    new_events: Dict[EventID, Event] = {}
    for dtag in reference_structure_dict:
        for event_id in event_dict:
            if event_id.dtag.dtag == dtag.dtag:
                new_events[event_id] = event_dict[event_id]

    return new_events


# ########
# Script
# #########

def main():
    # get args
    args: Args = Args.from_cmd()
    if args.debug > 0:
        Constants.DEBUG = args.debug

    # Get pandda directories
    system_path_list: List[Path] = list(path for path in args.pandda_dirs_dir.glob("*") if path.is_dir())
    if Constants.DEBUG > 0: print(f"Found {len(system_path_list)} systems")
    system_path_dict: SystemPathDict = SystemPathDict.from_dir(args.pandda_dirs_dir)

    # Get events tables: List[Path] -> Dict[EventID, Event]
    event_table_dict: Dict[System, pd.DataFrame] = get_event_table_dict(system_path_list)
    if Constants.DEBUG > 0: print(f"Found {len(event_table_dict)} system tables")

    # Get events
    event_dict: Dict[EventID, Event] = get_event_dict(event_table_dict, args.pandda_dirs_dir, args.autobuild_dirs_dir)
    if Constants.DEBUG > 0: print(f"Found {len(event_dict)} events")

    # Get reference models
    # reference_structure_dict: ReferenceStructureDict = ReferenceStructureDict.from_system_path_dict(system_path_dict)
    reference_structure_dict: ReferenceStructureDict = ReferenceStructureDict.from_structure_dir(
        args.reference_structure_dir)

    if Constants.DEBUG > 0: print(f"Found {len(reference_structure_dict)} reference models")

    # Filter events that have models
    event_with_reference_dict: Dict[EventID, Event] = get_event_with_reference_dict(
        event_dict,
        reference_structure_dict,
    )
    if Constants.DEBUG > 0: print(f"Found {len(event_with_reference_dict)} events with models")

    # Make output directory
    make_autobuild_output_dir(
        event_with_reference_dict,
        args.autobuild_dirs_dir,
    )

    # Map over all events
    map_dict(
        build_event,
        event_with_reference_dict,
    )


if __name__ == "__main__":
    main()
