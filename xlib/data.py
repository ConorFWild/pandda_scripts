from __future__ import annotations
from os import stat, system
import subprocess

from typing import *
from collections import abc
import dataclasses

import re
import json
from pathlib import Path

import numpy as np
import pandas as pd
import gemmi

from fragalysis_api.xcextracter.getdata import GetPdbData
from fragalysis_api.xcextracter.xcextracter import xcextracter


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

    PANDDA_JOB_LOG_FILE = "{system_name}.log"
    PANDDA_JOB_OUTPUT_FILE = "{system_name}.out"
    PANDDA_JOB_ERROR_FILE = "{system_name}.err"

    PANDDA_LOG_FILE = "pandda_log.json"
    
    ELBOW_COMMAND = "cd {event_dir.event_dir}; phenix. {smiles_file.smiles_file} --output=\"{autobuilding_ligand}\""
    FINAL_LIGAND_CIF_FILE = "final_ligand_cif.cif"
    
    MASKED_PDB_FILE = "masked.pdb"
    
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
    RHOFIT_SCRIPT_FILE = "run_rhofit.sh"
    
    PHASE_GRAFTED_MTZ_FILE = "phase_grafted_mtz.mtz"
    
    RHOFIT_HIT_REGEX = "(Hit_[^\s]+)[\s]+[^\s]+[\s]+[^\s]+[\s]+([^\s]+)"
    RHOFIT_CLUSTER_BUILD_REGEX = "Hit_([^_]+)_[^_]+_([^_]+).pdb"
    
    RHOFIT_EVENT_BUILD_RESULT_FILE = "result.json"
    
    CLUSTER = "CONDOR"
    CONDOR_JOB_ID_REGEX = r"[0-9]+\."
    CONDOR_STATUS_COMMAND = r"condor_q"
    
    XCHEM_DTAG_REGEX = r"^([0-9a-zA-Z]+[-]+[0-9a-zA-Z]+)"
    XCHEM_SYSTEM_REGEX = r"([^\-]+)\-[^\-]+"
    
    PHENIX_MAP_MODEL_CC_COMMAND = "phenix.map_model_cc {pdb_file} {event_map_file} resolution={resolution} ignore_symmetry_conflicts=True"
    PHENIX_MAP_MODEL_CC_LOG_FILE = "cc_per_residue.log"
    
    # PROJECT_CODE_MAPPING_DICT = {
    #     "NUDT7A"
    #     "ATAD"
    #     "BRD1A"
    #     "DCP2B"
    #     "FAM83BA"
    #     "MURD"
    #     "OXA100TA"
    #     "PARP14A"
    #     "PHIPA"
    #     "PTP1B"
    #     "SMTGR"
    #     "ATAD2A"
    #     "CAMK1DA"
    #     "DCLRE1AA"
    #     "FALZA"
    #     "HAO1A"
    #     "MUREECA"
    #     "NUDT21A"
    #     "NUDT4"
    #     "NUDT5A"
    #     "NUDT7A_CRUDE"
    #     "smTGRNEW"
    #     "STAG1A"
    #     "TBXTA"
    #     "VIM2"
    #     "XX02KALRNA"
    #     "TNCA"
    #     "ALAS2A"
    #     "EPB41L3A"
    #     "mArh"
    #     "INPP5DA"
    #     "nsp13"
    #     "Mac1"
    #     "Mpro"
    #     "NSP15_B"
    # }
    
    REFERENCE_STRUCTURE_FILE = "{dtag}.pdb"
    
    SKELETON_STRUCTURE_JSON_FILE="skeleton_score.json"
    SKELETON_XMAP_FILE_RECORD="xmap_file"
    SKELETON_STRUCTURE_FILE_RECORD="structure_file"
    SKELETON_SCORE_RECORD="score"
    
    
    

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
class ProjectCode:
    def __init__(self, project_code):
        self.project_code = project_code

    @staticmethod
    def from_pandda_processed_dir(pandda_dir: Path):
        processed_models_dir = pandda_dir / Constants.PANDDA_PROCESSED_DATASETS_DIR

        processed_model_dir = next(processed_models_dir.glob("*"))

        example_dtag = processed_model_dir.name

        project_name = System.from_dtag(example_dtag).system

        return ProjectCode(project_name)

    @staticmethod
    def from_dir(data_dir: Path):

        dtag_dir = next(data_dir.glob("*"))

        example_dtag = dtag_dir.name

        project_name = System.from_dtag(example_dtag).system

        return ProjectCode(project_name)


@dataclasses.dataclass()
class System:
    system: str
    
    def __hash__(self) -> int:
        return hash(self.system)
    
    @staticmethod
    def from_dtag(dtag: str):
        regex = Constants.XCHEM_SYSTEM_REGEX

        matches = re.findall(regex,
                             str(dtag),
                             )

        return System(matches[0])

@dataclasses.dataclass()
class Dtag:
    dtag: str
    
    def __hash__(self) -> int:
        return hash(self.dtag)
    
    @staticmethod
    def from_protein_code(protein_code):
        regex = Constants.XCHEM_DTAG_REGEX
        matches = re.findall(regex,
                             protein_code,
                             )

        return Dtag(matches[0])

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
        
    @staticmethod
    def from_system_row(system: System, row: pd.Series) -> EventID:
        dtag: Dtag = Dtag(row["dtag"])
        event_idx: EventIDX = EventIDX(row["event_idx"])
        return EventID(
            system=system,
            dtag=dtag,
            event_idx=event_idx,
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
    resolution: float
    
    @staticmethod
    def get_event_input_dir(system: System, dtag: Dtag, event_idx: EventIDX, pandda_dirs_dir: Path) -> Path:
        return pandda_dirs_dir / system.system /Constants.PANDDA_PROCESSED_DATASETS_DIR / dtag.dtag
    
    @staticmethod
    def get_event_output_dir(system: System, dtag: Dtag, event_idx: EventIDX, autobuild_dirs_dir: Path) -> Path:
        return autobuild_dirs_dir / system.system / dtag.dtag / event_idx.event_idx
    
    @staticmethod
    def from_system_row_datadir_autobuilddir(system: System, row: pd.Series, pandda_dirs_dir: Path, autobuild_dirs_dir: Path) -> Event:

        dtag: Dtag = Dtag(row["dtag"])
        event_idx: EventIDX = EventIDX(str(row["event_idx"]))
        bdc = row["1-BDC"]
        x = row["x"]
        y = row["y"]
        z = row["z"]
        resolution = row["analysed_resolution"]
        
        event_input_dir: Path = Event.get_event_input_dir(system, dtag, event_idx, pandda_dirs_dir)
        event_output_dir: Path = Event.get_event_output_dir(system, dtag, event_idx, autobuild_dirs_dir)
        
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
            resolution=resolution,
        )

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
class BuildID:
    system: System
    dtag: Dtag
    event_idx: EventIDX
    build_cluster: BuildClusterID
    build_number: BuildNumberID
    
    def __hash__(self) -> int:
        return hash(
            (
                self.system,
                self.dtag,
                self.event_idx,
                self.build_cluster,
                self.build_number,
                )
            )


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
            
        rhofit_dir:Path = event.event_output_dir / Constants.RHOFIT_DIR
        
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

@dataclasses.dataclass()
class EventTableDict:
    _dict: Dict[System, pd.DataFrame]
    
    def __getitem__(self, system: System) -> pd.DataFrame:
        return self._dict[system]
    
    def __setitem__(self, system: System, event_table: pd.DataFrame):
        self._dict[system] = event_table
        
    def __delitem__(self, system: System):
        del self._dict[system]
    
    def __iter__(self) -> Iterator[System]:
        for item in self._dict:
            yield item

    def __len__(self) -> int:
        return len(self._dict)
    
    @staticmethod
    def from_system_path_dict(system_path_dict: SystemPathDict) -> EventTableDict:
        event_table_dict = {}
        
        for system in system_path_dict:    
            # system: System = System(path.name)
            path: Path = system_path_dict[system]

            event_table_file: Path = path / Constants.PANDDA_ANALYSES_DIR / Constants.PANDDA_ANALYSE_EVENTS_FILE
            
            if Constants.DEBUG >0: print(event_table_file)
            
            if event_table_file.exists():
                event_table: pd.DataFrame = pd.read_csv(str(event_table_file))
                event_table_dict[system] = event_table
                
        return EventTableDict(event_table_dict)
        
    
    
@dataclasses.dataclass()
class EventDict:
    _dict: Dict[EventID, Event]
    
    def __getitem__(self, event_id: EventID) -> Event:
        return self._dict[event_id]
    
    def __setitem__(self, event_id: EventID, event: Event):
        self._dict[event_id] = event
        
    def __delitem__(self, event_id: EventID):
        del self._dict[event_id]
    
    def __iter__(self) -> Iterator[EventID]:
        for item in self._dict:
            yield item

    def __len__(self) -> int:
        return len(self._dict)
    
    @staticmethod
    def from_event_tables(event_table_dict: EventTableDict,
                          pandda_dirs_dir: Path, 
                          autobuild_dirs_dir: Path,
                          ) -> EventDict:
        event_dict: Dict[EventID, Event] = {}
    
        for system in event_table_dict:
            event_table: pd.DataFrame = event_table_dict[system]
            for index, row in event_table.iterrows():
                event_id: EventID = EventID.from_system_row(system, row)
                event: Event = Event.from_system_row_datadir_autobuilddir(
                    system, 
                    row, 
                    pandda_dirs_dir, 
                    autobuild_dirs_dir,
                    )
                
                event_dict[event_id] = event
                
        return EventDict(event_dict)
    
    def filter_reference_structure_dict(self, reference_structure_dict: ReferenceStructureDict) -> EventDict:
        new_events: Dict[EventID, Event] = {}
        
        for event_id in self:
            for dtag in reference_structure_dict:
                if dtag.dtag == event_id.dtag.dtag:
                    new_events[event_id] = self[event_id]

        return EventDict(new_events)    

@dataclasses.dataclass()
class Structure:
    structure: gemmi.Structure

    @staticmethod
    def from_model_path(path):
        structure = gemmi.read_structure(str(path))

        return Structure(structure)

    @staticmethod
    def from_pdb_string(pdb_string: str):
        structure = gemmi.read_pdb_string(pdb_string)
        return Structure(structure)
    
    def to_pdb(self, path: Path):
        self.structure.write_pdb(str(path))


@dataclasses.dataclass()
class Xmap:
    xmap: gemmi.FloatGrid

    @staticmethod
    def from_file(file):
        ccp4 = gemmi.read_ccp4_map(str(file))
        ccp4.setup()
        return Xmap(ccp4.grid)
    
    def to_grid(self) -> gemmi.FloatGrid:
        return self.xmap

@dataclasses.dataclass()
class SystemPathDict:
    _dict: Dict[System, Path]
    
    def __getitem__(self, key: System) -> Path:
        return self._dict[key]
    
    def __setitem__(self, key: System, value: Path):
        self._dict[key] = value
        
    def __delitem__(self, key: System):
        del self._dict[key]
    
    def __iter__(self) -> Iterator[System]:
        for item in self._dict:
            yield item

    def __len__(self) -> int:
        return len(self._dict)
    
    @staticmethod
    def from_dir(dir: Path) -> SystemPathDict:
        dir_list: List[Path] = list(path for path in dir.glob("*") if path.is_dir())
        system_dict_path: Dict[System, Path] = {System(dir.name): dir for dir in dir_list}
        
        return SystemPathDict(system_dict_path)


@dataclasses.dataclass()
class ReferenceStructureDict:
    _dict: Dict[Dtag, Structure]

    @staticmethod
    def from_system_path_dict(system_path_dict: SystemPathDict) -> ReferenceStructureDict:
        
        reference_structure_dict: Dict[Dtag, Structure] = {}
        for system in system_path_dict:
            system_path: Path = system_path_dict[system]
            print(f"Working on sytem: {system_path}")
            try:
                system_reference_structure_dict: Dict[Dtag, Structure] = ReferenceStructureDict.from_dir(system_path)
            except Exception as e:
                print(e)
                continue
            reference_structure_dict.update(system_reference_structure_dict)
            
        return ReferenceStructureDict(reference_structure_dict)
    
    def as_dict(self):
        return self._dict            
        

    @staticmethod
    def from_dir(pandda_dir: Path) -> Dict[Dtag, Structure]:
        print("\Getting project code...")
        project_code = ProjectCode.from_dir(pandda_dir)
        print("\tProject code is: {}".format(project_code.project_code))

        xcd = xcextracter(project_code.project_code)

        pdb_grabber = GetPdbData()

        structures = {}
        for index, row in xcd.iterrows():
            protein_code = row["protein_code"]
            dtag = Dtag.from_protein_code(protein_code)
            pdb_block = pdb_grabber.get_bound_pdb_file(protein_code)

            try:
                structure = Structure.from_pdb_string(pdb_block)
                structures[dtag] = structure

            except Exception as e:
                print(e)
                continue

        print("\tGot {} structures".format(len(structures)))
        return structures

    def to_dict(self) -> Dict[Dtag, Structure]:
        return self._dict

    def __iter__(self):
        for dtag in self._dict:
            yield dtag

    def __len__(self):
        return len(self._dict)

    def __getitem__(self, item):
        return self._dict[item]
    
    def __contains__(self, key) -> bool:
        return key in self._dict


@dataclasses.dataclass()
class BuildDict:
    _dict: Dict[BuildID, Build]
    
    def __getitem__(self, key: BuildID) -> Build:
        return self._dict[key]
    
    def __setitem__(self, key: BuildID, value: Build):
        self._dict[key] = value
        
    def __delitem__(self, key: BuildID):
        del self._dict[key]
    
    def __iter__(self) -> Iterator[BuildID]:
        for item in self._dict:
            yield item

    def __len__(self) -> int:
        return len(self._dict)
    
    @staticmethod
    def get_system_dict(autobuild_dir: Path):
        system_dict: Dict[System, Path] = {}
        
        system_path_list: List[Path] = list(autobuild_dir.glob("*"))
        
        for path in system_path_list:
            system: System = System(path.name)
            system_dict[system] = path

        return system_dict
    
    @staticmethod
    def get_dtag_path_dict(system_path_dict: Dict[System, Path]):
        dtag_path_dict: Dict[Tuple[System, Dtag], Path] = {}
        
        for system, path in system_path_dict.items():
            dtag_path_list: List[Path] = list(path.glob("*"))

            for dtag_path in dtag_path_list:
                dtag: Dtag = Dtag(dtag_path.name)
                dtag_path_dict[(system, dtag)] = dtag_path

        return dtag_path_dict

    @staticmethod
    def get_event_path_dict(dataset_path_dict: Dict[Tuple[System, Dtag], Path]):
        event_path_dict: Dict[Tuple[System, Dtag, EventIDX], Path] = {}
        
        for system_dtag, path in dataset_path_dict.items():
            system: System = system_dtag[0]
            dtag: Dtag = system_dtag[1]
            event_path_list: List[Path] = list(path.glob("*"))

            for event_path in event_path_list:
                event_idx: EventIDX = EventIDX(event_path.name)
                event_path_dict[(system, dtag, event_idx)] = event_path

        return event_path_dict
    
    @staticmethod
    def get_build_dict(dataset_path_dict: Dict[Tuple[System, Dtag, EventIDX], Path]) -> Dict[BuildID, Build]:
        build_dict: Dict[BuildID, Build] = {}
        
        for system_dtag_eventidx, path in dataset_path_dict.items():
            system: System = system_dtag_eventidx[0]
            dtag: Dtag = system_dtag_eventidx[1]
            event_idx: EventIDX = system_dtag_eventidx[2]

            build_result_file: Path = path / Constants.RHOFIT_RESULT_JSON_FILE
            
            if not build_result_file.exists():
                continue
            
            with open(build_result_file, "r") as f:
                result_dict = json.load(f)

            # event_build_results: EventBuildResults = EventBuildResults.from_dict(result_dict)
            
            for build_cluster_id_str in result_dict:
                build_cluster_id = BuildClusterID(build_cluster_id_str)
                
                build_cluster_dict: Dict[str, Dict[str, str]] = result_dict[build_cluster_id_str]
                
                for build_number_id_str in build_cluster_dict:
                    current_build_dict: Dict[str, str] = build_cluster_dict[build_number_id_str]
                    build_file: Path = Path(current_build_dict["build_file"])
                    build_rscc: float = float(current_build_dict["build_rscc"])
                    build: Build = Build(build_file, build_rscc)
                    build_number_id: BuildNumberID = BuildNumberID(int(build_number_id_str))
                    build_id: BuildID = BuildID(system, dtag, event_idx, build_cluster_id, build_number_id)
                    build_dict[build_id] = build

        return build_dict
     
                
    
    @staticmethod
    def from_autobuild_dir(autobuild_dir: Path) -> BuildDict:
        
        system_path_dict: Dict[System, Path] = BuildDict.get_system_dict(autobuild_dir)
        dataset_path_dict: Dict[Tuple[System, Dtag], Path] = BuildDict.get_dtag_path_dict(system_path_dict)
        event_path_dict : Dict[Tuple[System, Dtag, EventIDX], Path] = BuildDict.get_event_path_dict(dataset_path_dict)
        build_number_dict : Dict[BuildID, Build] = BuildDict.get_build_dict(event_path_dict)
        
        return BuildDict(build_number_dict)

@dataclasses.dataclass()
class LigandResidues:
    residues: List[Any]

    def __iter__(self):
        for residue in self.residues:
            yield residue

    @staticmethod
    def from_structure(structure: Structure):
        residues = []
        for model in structure.structure:
            for chain in model:
                for residue in chain:
                    if residue.name == "LIG":
                        residues.append(residue)

        return LigandResidues(residues)

@dataclasses.dataclass()
class RMSD:
    rmsd: float

    @staticmethod
    def from_structures(reference_structure: Structure, build_structure: Structure) -> RMSD:
        
        reference_ligand = LigandResidues.from_structure(reference_structure)
        autobuild_ligand = LigandResidues.from_structure(build_structure)

        residue_rmsds = []
        for residue_1 in autobuild_ligand.residues:
            for residue_2 in reference_ligand.residues:
                rmsd = RMSD.from_residues(residue_1,
                                            residue_2,
                                            )

                residue_rmsds.append(rmsd.rmsd)

        rmsd = RMSD(min(residue_rmsds))
            
        return rmsd

    @staticmethod
    def from_residues(residue_1, residue_2):
        closest_distances = []
        for atom_1 in residue_1:
            distances = []
            for atom_2 in residue_2:
                distance = atom_1.pos.dist(atom_2.pos)
                distances.append(distance)

            closest_distance = min(distances)
            closest_distances.append(closest_distance)

        mean_distance = float(np.mean(closest_distances))
        return RMSD(mean_distance)

@dataclasses.dataclass()
class RMSDDict:
    _dict: Dict[BuildID, RMSD]
    
    
    def __getitem__(self, key: BuildID) -> RMSD:
        return self._dict[key]
    
    def __setitem__(self, key: BuildID, value: RMSD):
        self._dict[key] = value
        
    def __delitem__(self, key: BuildID):
        del self._dict[key]
    
    def __iter__(self) -> Iterator[BuildID]:
        for item in self._dict:
            yield item

    def __len__(self) -> int:
        return len(self._dict)
    
    @staticmethod
    def from_build_dict(
        build_dict: BuildDict,
        reference_structure_dict: ReferenceStructureDict,
        ) -> RMSDDict:
        rmsd_dict: Dict[BuildID, RMSD] = {}
        for build_id in build_dict:
            build: Build = build_dict[build_id]
            dtag: Dtag = build_id.dtag
            
            if dtag in reference_structure_dict:
                reference_structure: Structure = reference_structure_dict[dtag]
                
                build_structure_file: Path = build.build_file
                
                build_structure: Structure = Structure.from_model_path(build_structure_file)
                rmsd: RMSD = RMSD.from_structures(reference_structure, build_structure)
                
                rmsd_dict[build_id] = rmsd

        return RMSDDict(rmsd_dict)            
    
    def best_by_dtag(self) -> RMSDDict:
        best_rmsd_dict: Dict[BuildID, RMSD] = {}
        
        dtag_to_buildid_dict: Dict[Dtag, List[BuildID]] = {}
        for build_id in self:
            dtag: Dtag = build_id.dtag
            if dtag not in dtag_to_buildid_dict:
                dtag_to_buildid_dict[dtag] = []
            dtag_to_buildid_dict[dtag].append(build_id)
        
        for dtag in dtag_to_buildid_dict:
            build_id_list: List[BuildID] = dtag_to_buildid_dict[dtag]
            rmsd_list: List[RMSD] = [self[build_id] for build_id in build_id_list]
            rmsd_float_list: List[float] = [rmsd.rmsd for rmsd in rmsd_list]
            
            min_rmsd_index: int = np.argmin(rmsd_float_list)
            
            best_build_id: BuildID = build_id_list[min_rmsd_index]
            best_rmsd: RMSD = rmsd_list[min_rmsd_index]

            best_rmsd_dict[best_build_id] = best_rmsd
            
        return RMSDDict(best_rmsd_dict)
            
def get_residue_centroid(residue: gemmi.Residue) -> Tuple[float, float, float]:
    pos_list = []
    for atom in residue:
        pos = atom.pos
        pos_tuple = (pos.x, pos.y, pos.z)
        pos_list.append(pos_tuple)
        
    pos_array = np.array(pos_list)
    centroid = np.mean(pos_array, axis=0)
    centroid_tuple = (centroid[0], centroid[1], centroid[2])
    return centroid_tuple
            
def get_event_distance_from_reference_model_dict(
    event_dict: EventDict, 
    reference_structure_dict: ReferenceStructureDict,
    system_path_dict: SystemPathDict,
    ) -> Dict[Dtag, float]:
    dtag_distance_from_reference_model_dict: Dict[Dtag, float] = {}
    
    # Iterate over reference models
    for dtag in reference_structure_dict:
        # Get system
        system: System = System.from_dtag(dtag.dtag)
        # Check if the pandda finished
        csv_path: Path = system_path_dict[system] / Constants.PANDDA_ANALYSES_DIR / Constants.PANDDA_ANALYSE_EVENTS_FILE
        # Continue if not finished
        if not csv_path.exists():
            print(f"System didn't finish: {system}")
            continue
        # Get dtag path
        dtag_path: Path = system_path_dict[system] / Constants.PANDDA_PROCESSED_DATASETS_DIR / dtag.dtag
        # Check if this system has been processed, and skip if not
        if not dtag_path.exists():
            print(f"Dtag has not been processed: {dtag_path}, skipping")
            continue
        # Get reference structure
        reference_structure: Structure = reference_structure_dict[dtag]
        # Get ligands
        ligand_residues: LigandResidues = LigandResidues.from_structure(reference_structure)
        # Iterate over events with this dtag
        dtag_event_list: List[EventID] = [event_id for event_id in event_dict if event_id.dtag == dtag]

        if len(dtag_event_list) == 0:
            dtag_distance_from_reference_model_dict[dtag] = -1
            continue


        # Iterate
        ligand_distance_list: List[float] = []
        for ligand_residue in ligand_residues:
            # Get ligand centroid
            ligand_centroid: Tuple[float, float, float] = get_residue_centroid(ligand_residue)                
            
            for event_id in dtag_event_list:
                # Get event
                event: Event = event_dict[event_id]
                # Get event centroid
                event_centroid: Tuple[float, float, float] = (event.x, event.y, event.z)
                # Get distance
                distance: float = np.linalg.norm(np.array(event_centroid) - np.array(ligand_centroid))
                # Update
                ligand_distance_list.append(distance)
                
        # Skip processing if none found
        if len(ligand_distance_list) == 0:
            continue
        # Get minimum
        min_distance = min(ligand_distance_list)
        # update
        dtag_distance_from_reference_model_dict[dtag] = min_distance
        
    return dtag_distance_from_reference_model_dict
                

    
# @dataclasses.dataclass()
# class RSCCDict:
#     _dict: Dict[EventID, RSCC]
    
#     @staticmethod
#     def from_build_dict() -> RSCCDict:
        
@dataclasses.dataclass()
class RSCCDict:
    _dict: Dict[EventID, float]
    
    def __getitem__(self, key: BuildID) -> RMSD:
        return self._dict[key]
    
    def __setitem__(self, key: BuildID, value: RMSD):
        self._dict[key] = value
        
    def __delitem__(self, key: BuildID):
        del self._dict[key]
    
    def __iter__(self) -> Iterator[BuildID]:
        for item in self._dict:
            yield item

    def __len__(self) -> int:
        return len(self._dict)
    
    @staticmethod
    def get_phenix_map_model_cc_log_file(reference_structure_file,
            event_map_file,
            resolution,
            event_analysis_dir: Path,
            ) -> Path:
        # Path to ligand cc log
        ligand_cc_log_file: Path = event_analysis_dir / Constants.PHENIX_MAP_MODEL_CC_LOG_FILE
        
        # format script
        command: str = Constants.PHENIX_MAP_MODEL_CC_COMMAND.format(
            dir=str(event_analysis_dir),
            pdb_file=str(reference_structure_file),
            event_map_file=str(event_map_file),
            resolution=resolution,
        )
        
        # Run
        p = subprocess.Popen(
            command,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.subprocess.PIPE,
        )
        
        stdout, stderr = p.communicate()
        
        return ligand_cc_log_file
    
    @staticmethod
    def get_phenix_map_model_cc_log_string(phenix_map_model_cc_log_file: Path) -> str:
        with open(str(phenix_map_model_cc_log_file), "r") as f:
            string = f.read()
            
        return string
    
    @staticmethod
    def parse_phenix_map_model_cc_log(phenix_map_model_cc_log_string: str) -> float:
        pattern: Pattern = Constants.PHENIX_MAP_MODEL_CC_pattern
        matches = re.findall(
            pattern,
            phenix_map_model_cc_log_string,
            )
        return float(matches[0])
            
    @staticmethod
    def get_rscc(
        event: Event,
        reference_structure: Structure,
        pandda_dirs_dict: SystemPathDict,
        analyses_dir: Path,
        ):
        # Get analysis dir
        event_analysis_dir: Path = analyses_dir / event.system.system / event.dtag.dtag / event.event_idx.event_idx
        
        # Write reference structure
        reference_structure_file: Path = reference_structure.to_pdb(
            event_analysis_dir,
            )
        
        # Get PanDDA path
        for system in pandda_dirs_dir:
            if system.system == event.system.system:
                pandda_dir: Path = pandda_dirs_dict[system]
        
        
        # Get event map file
        event_map_file: Path = pandda_dir / Constants.PANDDA_PROCESSED_DATASETS_DIR / event.dtag.dtag / Constants.PANDDA_EVENT_MAP_FILE.format(
            dtag=event.dtag.dtag, 
            event_idx=event.event_idx.event_idx,
            bdc=event.bdc,
            )
        
        # Run RSCC calc
        phenix_map_model_cc_log_file: Path = RSCCDict.get_phenix_map_model_cc_log_file(
            reference_structure_file,
            event_map_file,
            event.resolution,
            event_analysis_dir,
        )
        
        # Read log
        phenix_map_model_cc_log_string: str = RSCCDict.get_phenix_map_model_cc_log_string(
            event_map_file
        )
        
        # Parse log
        rscc: float = RSCCDict.parse_phenix_map_model_cc_log(phenix_map_model_cc_log_string)
        
        return rscc
        
        
    
    @staticmethod
    def from_event_dict(
        event_dict: EventDict,
        reference_structure_dict: ReferenceStructureDict, 
        pandda_dirs_dir: SystemPathDict,
        analyses_dir: Path,
        ) -> RSCCDict:
        
        event_rscc_dict: Dict[EventID, float] = {}
        
        for event_id in event_dict:
            # Check if there is a corresponding reference structure for the event, and get rscc if so
            for dtag in reference_structure_dict:
                if dtag.dtag == event_id.dtag.dtag:
                    rscc: float = RSCCDict.get_rscc(
                        event_dict[event_id],
                        reference_structure_dict[dtag],
                        pandda_dirs_dir,
                        analyses_dir,
                        )
                    event_dict[event_id] = rscc
        return RSCCDict(event_rscc_dict)
    
    def filter_best(self) -> RSCCDict:
        dtag_rscc_dict: Dict[Dtag, List[float]] = {}
        for event_id in self:
            if event_id.dtag not in dtag_rscc_dict:
                dtag_rscc_dict[event_id.dtag] = []
            
            dtag_rscc_dict[event_id.dtag].append(self[event_id])
        
        best_rscc_dict: Dict[Dtag, float] = {}
        # Get min
        for dtag in dtag_rscc_dict:
            best_rscc_dict[dtag] = min(dtag_rscc_dict[dtag])
            
        return RSCCDict(best_rscc_dict)
    

@dataclasses.dataclass()
class Node:
    x: float
    y: float
    z: float

@dataclasses.dataclass()
class Skeleton:
    _node_list: List[Node] 
    
    def __iter__(self):
        for node in self._node_list:
            yield node
    
    @staticmethod
    def from_structure(_structure: Structure) -> Skeleton:
        
        structure = _structure.structure
        
        # Get ligand res
        for n_ch, chain in enumerate(structure[0]):
            for n_res, res in enumerate(chain):
                if res.name == "LIG":
                    residue = res
                    continue
        
        # Make the neighbourhood search
        ns = gemmi.NeighborSearch(structure[0], structure.cell, 3)
        for n_ch, chain in enumerate(structure[0]):
            for n_res, res in enumerate(chain):
                for n_atom, atom in enumerate(res):
                    if not atom.is_hydrogen():
                        ns.add_atom(atom, n_ch, n_res, n_atom)

        # Find the bonds
        bond_list = []
        numbered_atoms = {}
        for atom_idx, atom in enumerate(residue):
            numbered_atoms[atom_idx] = atom
            # A little more than the average distance of a standard c-c bond
            marks = ns.find_neighbors(atom, min_dist=0.1, max_dist=1.6)
            for mark in marks:
                key = (atom_idx, mark.atom_idx)
                bond_list.append(key)
        print(f"Found {len(numbered_atoms)} residue atoms")
                
        # Filter duplicates
        filtered_bonds = []
        for bond in bond_list:
            reversed_bond = (bond[1], bond[0])
            if reversed_bond not in filtered_bonds:
                filtered_bonds.append(bond)
        print(f"Found {len(filtered_bonds)} bonds")

        # Get normal atom
        node_list = []
        for atom in residue:
            pos = atom.pos
            node = Node(
                pos.x,
                pos.y,
                pos.z,
            )
            node_list.append(node)
        
        # Get bonds
        for bond in filtered_bonds:
            atom_1 = numbered_atoms[bond[0]]
            atom_2 = numbered_atoms[bond[1]]
            
            pos_1 = atom_1.pos
            pos_2 = atom_2.pos
            
            mean_x = (pos_1.x + pos_2.x) / 2        
            mean_y = (pos_1.y + pos_2.y) / 2        
            mean_z = (pos_1.z + pos_2.z) / 2        

            node = Node(
                mean_x,
                mean_y,
                mean_z,
            )
            node_list.append(node)
            
        print(f"Found {len(node_list)} nodes")
            
        return Skeleton(node_list)
    
    
@dataclasses.dataclass()
class SkeletonScore:
    skeleton_score: float
    
    @staticmethod
    def from_structure(structure: Structure,
                     xmap: Xmap,
                     contour: float=1.0,
                     ) -> SkeletonScore:
        
        
        skeleton: Skeleton = Skeleton.from_structure(structure)
        
        xmap_grid: gemmi.FloatGrid = xmap.to_grid()
        
        xmap_array: np.array = np.array(xmap_grid, copy=False)
        
        xmap_array[xmap_array < contour] = 0
        
        values: List[float] = []
        for node in skeleton:
            pos = gemmi.Position(node.x,
                                 node.y,
                                 node.z,
                                 )
            
            sample_value: float = xmap_grid.interpolate_value(pos)
            
            if sample_value > 0:
                values.append(1)
            else:
                values.append(0)
                
        return SkeletonScore(sum(values) / len(values))

