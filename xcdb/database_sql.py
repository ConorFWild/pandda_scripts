import argparse
import dataclasses
import json
import os
import subprocess
from pathlib import Path
from typing import *

import gemmi
import numpy as np
import pandas as pd
from sqlalchemy import (Boolean, Column, Float, ForeignKey, Integer, String,
                        create_engine, func)
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship, sessionmaker


import fire


class Constants:
    COMPOUND_TABLE = "compound"
    SMILES_TABLE = "smiles"
    REFLECTIONS_TABLE = "reflections"
    RESOLUTION_TABLE = "resolution"
    SPACEGROUP_TABLE = "spacegroup"
    UNIT_CELL_TABLE = "unit_cell"
    MODEL_TABLE = "model"
    DATASET_TABLE = "dataset"
    SYSTEM_TABLE = "system"
    SEQUENCE_IDENTITY_TABLE = "system_identity"
    PANDDA_TABLE = "pandda"
    PANDDA_ERROR_TABLE = "pandda_error"
    EVENT_TABLE = "event"
    REFERENCE_TABLE = "reference"
    AUTOBUILD_TABLE = "autobuild"
    AUTOBUILD_BEST_TABLE = "autobuild_best"
    AUTOBUILD_SKELETON_SCORE_TABLE = "autobuild_skeleton_score"
    AUTOBUILD_RSCC_TABLE = "autobuild_rscc"
    AUTOBUILD_RMSD_TABLE = "autobuild_rmsd"
    AUTOBUILD_HUMAN_TABLE = "autobuild_human"
    REAL_SPACE_CLUSTERING_TABLE = "real_space_clustering"
    EVENT_SCORE_TABLE = "event_score"
    TEST_TABLE = "test"
    TRAIN_TABLE = "train"


base = declarative_base()


class Compound(base):
    __tablename__ = Constants.COMPOUND_TABLE
    id = Column(Integer, primary_key=True)
    path = Column(String(255))


class Smiles(base):
    __tablename__ = Constants.SMILES_TABLE
    id = Column(Integer, primary_key=True)
    path = Column(String(255))


class Reflections(base):
    __tablename__ = Constants.REFLECTIONS_TABLE
    id = Column(Integer, primary_key=True)
    path = Column(String(255))


class Resolution(base):
    __tablename__ = Constants.RESOLUTION_TABLE
    id = Column(Integer, primary_key=True)
    resolution = Column(Float)

    reflections_id = Column(Integer, ForeignKey(Reflections.id))
    reflections = relationship(Reflections)


class Spacegroup(base):
    __tablename__ = Constants.SPACEGROUP_TABLE
    id = Column(Integer, primary_key=True)
    spacegroup = Column(Integer)

    reflections_id = Column(Integer, ForeignKey(Reflections.id))
    reflections = relationship(Reflections)


class UnitCell(base):
    __tablename__ = Constants.UNIT_CELL_TABLE
    id = Column(Integer, primary_key=True)
    a = Column(Float)
    b = Column(Float)
    c = Column(Float)
    alpha = Column(Float)
    beta = Column(Float)
    gamma = Column(Float)

    reflections_id = Column(Integer, ForeignKey(Reflections.id))
    reflections = relationship(Reflections)


class Model(base):
    __tablename__ = Constants.MODEL_TABLE
    id = Column(Integer, primary_key=True)
    path = Column(String(255))


class System(base):
    __tablename__ = Constants.SYSTEM_TABLE
    id = Column(Integer, primary_key=True)
    system = Column(String(255))
    path = Column(String(255))


class SequenceIdentity(base):
    __tablename__ = Constants.SEQUENCE_IDENTITY_TABLE
    id = Column(Integer, primary_key=True)
    identity = Column(Float)

    system_1_id = Column(Integer, ForeignKey(System.id))
    system_2_id = Column(Integer, ForeignKey(System.id))

    system_1 = relationship(System, foreign_keys=[system_1_id])
    system_2 = relationship(System, foreign_keys=[system_2_id])


class Dataset(base):
    __tablename__ = Constants.DATASET_TABLE
    id = Column(Integer, primary_key=True)
    dtag = Column(String(255))
    path = Column(String(255))

    # Foriegn keys
    system_id = Column(Integer, ForeignKey(System.id))
    reflections_id = Column(Integer, ForeignKey(Reflections.id))
    model_id = Column(Integer, ForeignKey(Model.id))
    compound_id = Column(Integer, ForeignKey(Compound.id))
    smiles_id = Column(Integer, ForeignKey(Smiles.id))

    # Relationships
    system = relationship(System)
    reflections = relationship(Reflections)
    model = relationship(Model)
    compound = relationship(Compound)
    smiles = relationship(Smiles)


class PanDDA(base):
    __tablename__ = Constants.PANDDA_TABLE
    id = Column(Integer, primary_key=True)
    success = Column(Boolean)
    runtime = Column(Float)
    path = Column(String(255))

    # Foreign keys
    system_id = Column(Integer, ForeignKey(System.id))

    # Relationships
    system = relationship(System)


class PanDDAError(base):
    __tablename__ = Constants.PANDDA_ERROR_TABLE
    id = Column(Integer, primary_key=True)
    path = Column(String(255))
    error = Column(String(2047))

    # Foreign keys
    system_id = Column(Integer, ForeignKey(System.id))
    pandda_id = Column(Integer, ForeignKey(PanDDA.id))

    # Relationships
    system = relationship(System)
    pandda = relationship(PanDDA)


class Event(base):
    __tablename__ = Constants.EVENT_TABLE
    id = Column(Integer, primary_key=True)
    event_idx = Column(Integer)
    x = Column(Float)
    y = Column(Float)
    z = Column(Float)
    analysed_Resolution = Column(Float)
    bdc = Column(Float)
    event_map = Column(String(255))

    # Foreign keys
    system_id = Column(Integer, ForeignKey(System.id))
    dataset_id = Column(Integer, ForeignKey(Dataset.id))
    pandda_id = Column(Integer, ForeignKey(PanDDA.id))

    # Relationships
    system = relationship(System)
    dataset = relationship(Dataset)
    pandda = relationship(PanDDA)


class ReferenceModel(base):
    __tablename__ = Constants.REFERENCE_TABLE
    id = Column(Integer, primary_key=True)
    path = Column(String(255))

    # Foreign keys
    system_id = Column(Integer, ForeignKey(System.id))
    dataset_id = Column(Integer, ForeignKey(Dataset.id))

    # Relationships
    system = relationship(System)
    dataset = relationship(Dataset)


class Autobuild(base):
    __tablename__ = Constants.AUTOBUILD_TABLE
    id = Column(Integer, primary_key=True)
    path = Column(String(255))

    # Foreign keys
    system_id = Column(Integer, ForeignKey(System.id))
    dataset_id = Column(Integer, ForeignKey(Dataset.id))
    pandda_id = Column(Integer, ForeignKey(PanDDA.id))
    event_id = Column(Integer, ForeignKey(Event.id))

    # Relationships
    system = relationship(System)
    dataset = relationship(Dataset)
    pandda = relationship(PanDDA)
    event = relationship(Event)


class AutobuildBest(base):
    __tablename__ = Constants.AUTOBUILD_BEST_TABLE
    id = Column(Integer, primary_key=True)
    path = Column(String(255))

    # Foreign keys
    system_id = Column(Integer, ForeignKey(System.id))
    dataset_id = Column(Integer, ForeignKey(Dataset.id))
    pandda_id = Column(Integer, ForeignKey(PanDDA.id))
    event_id = Column(Integer, ForeignKey(Event.id))
    autobuild_id = Column(Integer, ForeignKey(Autobuild.id))

    # Relationships
    system = relationship(System)
    dataset = relationship(Dataset)
    pandda = relationship(PanDDA)
    event = relationship(Event)
    autobuild = relationship(Autobuild)


class AutobuildSkeletonScore(base):
    __tablename__ = Constants.AUTOBUILD_SKELETON_SCORE_TABLE
    id = Column(Integer, primary_key=True)
    autobuild_id = Column(Integer, ForeignKey(Autobuild.id))
    skeleton_score = Column(Float)

    # Foreign keys
    system_id = Column(Integer, ForeignKey(System.id))
    dataset_id = Column(Integer, ForeignKey(Dataset.id))
    pandda_id = Column(Integer, ForeignKey(PanDDA.id))
    event_id = Column(Integer, ForeignKey(Event.id))
    autobuild_id = Column(Integer, ForeignKey(Autobuild.id))

    # Relationships
    system = relationship(System)
    dataset = relationship(Dataset)
    pandda = relationship(PanDDA)
    event = relationship(Event)
    autobuild = relationship(Autobuild)


class AutobuildRSCC(base):
    __tablename__ = Constants.AUTOBUILD_RSCC_TABLE
    id = Column(Integer, primary_key=True)
    rscc = Column(Float)

    # Foreign keys
    system_id = Column(Integer, ForeignKey(System.id))
    dataset_id = Column(Integer, ForeignKey(Dataset.id))
    pandda_id = Column(Integer, ForeignKey(PanDDA.id))
    event_id = Column(Integer, ForeignKey(Event.id))
    autobuild_id = Column(Integer, ForeignKey(Autobuild.id))

    # Relationships
    system = relationship(System)
    dataset = relationship(Dataset)
    pandda = relationship(PanDDA)
    event = relationship(Event)
    autobuild = relationship(Autobuild)


class AutobuildRMSD(base):
    __tablename__ = Constants.AUTOBUILD_RMSD_TABLE
    id = Column(Integer, primary_key=True)
    reference_id = Column(Integer, ForeignKey(ReferenceModel.id))
    rmsd = Column(Float)

    # Foreign keys
    system_id = Column(Integer, ForeignKey(System.id))
    dataset_id = Column(Integer, ForeignKey(Dataset.id))
    pandda_id = Column(Integer, ForeignKey(PanDDA.id))
    event_id = Column(Integer, ForeignKey(Event.id))
    autobuild_id = Column(Integer, ForeignKey(Autobuild.id))

    # Relationships
    system = relationship(System)
    dataset = relationship(Dataset)
    pandda = relationship(PanDDA)
    event = relationship(Event)
    autobuild = relationship(Autobuild)


class HumanAutobuildComparison(base):
    __tablename__ = Constants.AUTOBUILD_HUMAN_TABLE
    id = Column(Integer, primary_key=True)

    # foreign keys
    system_id = Column(Integer, ForeignKey(System.id))
    dataset_id = Column(Integer, ForeignKey(Dataset.id))
    pandda_id = Column(Integer, ForeignKey(PanDDA.id))
    event_id = Column(Integer, ForeignKey(Event.id))
    autobuild_id = Column(Integer, ForeignKey(Autobuild.id))
    reference_id = Column(Integer, ForeignKey(ReferenceModel.id))

    # Relationships
    system = relationship(System)
    dataset = relationship(Dataset)
    pandda = relationship(PanDDA)
    event = relationship(Event)
    autobuild = relationship(Autobuild)
    reference_model = relationship(ReferenceModel)


class RealSpaceClustering(base):
    __tablename__ = Constants.REAL_SPACE_CLUSTERING_TABLE
    id = Column(Integer, primary_key=True)
    table_path = Column(String(255))
    graph_path = Column(String(255))

    # Foreign keys
    system_id = Column(Integer, ForeignKey(System.id))

    # Relationships
    system = relationship(System)


class EventScore(base):
    __tablename__ = Constants.EVENT_SCORE_TABLE
    id = Column(Integer, primary_key=True)
    event_score = Column(Float)

    # Foreign keys
    system_id = Column(Integer, ForeignKey(System.id))
    dataset_id = Column(Integer, ForeignKey(Dataset.id))
    pandda_id = Column(Integer, ForeignKey(PanDDA.id))
    event_id = Column(Integer, ForeignKey(Event.id))

    # Relationships
    system = relationship(System)
    dataset = relationship(Dataset)
    pandda = relationship(PanDDA)
    event = relationship(Event)


class Test(base):
    __tablename__ = Constants.TEST_TABLE
    id = Column(Integer, primary_key=True)

    system_id = Column(Integer, ForeignKey(System.id))
    system = relationship(System)


class Train(base):
    __tablename__ = Constants.TRAIN_TABLE
    id = Column(Integer, primary_key=True)

    system_id = Column(Integer, ForeignKey(System.id))
    system = relationship(System)


class Database:

    def __init__(self, database_path: Path, overwrite: bool = False) -> None:
        super().__init__()
        # conn = sqlite3.connect('example.db')
        if overwrite:
            if database_path.exists():
                os.remove(str(database_path))

        engine = create_engine(f'sqlite:///{str(database_path)}')

        base.metadata.bind = engine
        base.metadata.create_all(engine)

        DBSession = sessionmaker()
        DBSession.bind = engine
        session = DBSession()

        self.session = session

    def populate_systems(self, system_dir: Path):
        system_path_list: List[Path] = list(path for path in system_dir.glob("*") if path.is_dir())

        for system_path in system_path_list:
            system_name: str = system_path.name
            system_path: str = str(system_path)

            system = System(
                system=system_name,
                path=system_path,
            )
            self.session.add(system)

        self.session.commit()

        print("Populated systems")
        print(f"Got {self.session.query(func.count(System.id)).scalar()} systems")

    def populate_sequence_identity(self, sequence_identity_dir: Path):
        sequence_identity_file = sequence_identity_dir / data.Constants.IDENTITY_TABLE_FILE

        table = pd.read_csv(str(sequence_identity_file))

        for index, row in table.iterrows():
            system_name_1 = row["system"]

            for column_index, value in row.items():
                if column_index == "system":
                    continue

                system_name_2 = column_index

                identity = value

                system_1 = self.session.query(System).filter(System.system == system_name_1).first()
                system_2 = self.session.query(System).filter(System.system == system_name_2).first()

                sequence_identity = SequenceIdentity(identity=identity,
                                                     system_1=system_1,
                                                     system_2=system_2,
                                                     )

                self.session.add(sequence_identity)

        self.session.commit()

        print("Populated sequence identity")
        print(f"Got {self.session.query(func.count(SequenceIdentity.id)).scalar()} systems")

    def populate_models_reflections_compounds_smiles_datasets(self):

        for system in self.session.query(System).all():
            system_path = Path(system.path)

            # dataset_dirs_dir = system_path / system.system
            # print(dataset_dirs_dir)

            dataset_dir_list = list(path for path in system_path.glob("*") if path.is_dir())

            print(f"Got {len(dataset_dir_list)} data dirs for system {system.system}")

            for dataset_dir in dataset_dir_list:
                # compound_path = dataset_dir / Constants.COMPOUND_FILE
                # reflections_path = dataset_dir / Constants.REFLECTIONS_FILE
                # model_path = dataset_dir / Constants.MODEL_FILE

                dataset_compounds = list(dataset_dir.glob("*.cif"))
                if len(dataset_compounds) == 0:
                    compound_path = None
                else:
                    compound_path = dataset_compounds[0]

                dataset_reflections = list(dataset_dir.glob("*.mtz"))
                if len(dataset_reflections) == 0:
                    reflections_path = None
                else:
                    reflections_path = dataset_reflections[0]

                dataset_models = list(dataset_dir.glob("*.pdb"))
                if len(dataset_models) == 0:
                    model_path = None
                else:
                    model_path = dataset_models[0]

                dataset_smiles = list(dataset_dir.glob("*.smiles"))
                if len(dataset_smiles) == 0:
                    smiles_path = None
                else:
                    smiles_path = dataset_smiles[0]

                model = Model(path=str(model_path))

                reflections = Reflections(path=str(reflections_path))

                compound = Compound(path=str(compound_path))

                smiles = Smiles(path=str(smiles_path))

                # system_name = xlib.data.System.from_dtag(dataset_dir.name)   
                # print(system_name)         

                # system_obj = self.session.query(System).filter(System.system == system_name.system).first()

                dataset = Dataset(dtag=dataset_dir.name,
                                  path=str(dataset_dir),
                                  system=system,
                                  reflections=reflections,
                                  compound=compound,
                                  smiles=smiles,
                                  model=model,
                                  )

                self.session.add(model)
                self.session.add(reflections)
                self.session.add(compound)
                self.session.add(smiles)
                self.session.add(dataset)

        self.session.commit()

        print("Populated systems")
        print(
            (
                f"Got {self.session.query(func.count(Model.id)).scalar()} models\n"
                f"Got {self.session.query(func.count(Reflections.id)).scalar()} Reflections\n"
                f"Got {self.session.query(func.count(Compound.id)).scalar()} Compounds\n"
                f"Got {self.session.query(func.count(Dataset.id)).scalar()} datasets\n"
            )
        )

    def populate_reference_models(self, reference_structure_dir: Path):

        reference_model_list = reference_structure_dir.glob("*")

        for reference_model_path in reference_model_list:

            dtag = reference_model_path.stem

            # system = xlib.data.System.from_dtag(dtag)            

            # system = self.session.query(System).filter(System.system == system.system).first()

            dataset = self.session.query(Dataset).filter(Dataset.dtag.contains(dtag)).first()

            if not dataset:
                continue

            system = dataset.system

            reference_model = ReferenceModel(
                path=str(reference_model_path),
                system=system,
                dataset=dataset,
            )
            self.session.add(reference_model)

        self.session.commit()

        print("Populated reference models")
        print(
            (
                f"Got {self.session.query(func.count(ReferenceModel.id)).scalar()} models\n"
            )
        )

    def populate_resolution_spacegroup_unit_cell(self):
        for reflections in self.session.query(Reflections).all():
            try:
                reflections_path = Path(reflections.path)

                mtz = gemmi.read_mtz_file(str(reflections_path))

                res = mtz.resolution_high()

                sg = mtz.spacegroup.ccp4

                resolution = Resolution(
                    resolution=res,
                    reflections=reflections,
                )
                self.session.add(resolution)

                spacegroup = Spacegroup(
                    spacegroup=sg,
                    reflections=reflections,
                )
                self.session.add(spacegroup)

                unit_cell = UnitCell(
                    a=mtz.cell.a,
                    b=mtz.cell.b,
                    c=mtz.cell.c,
                    alpha=mtz.cell.alpha,
                    beta=mtz.cell.beta,
                    gamma=mtz.cell.gamma,
                )
                self.session.add(unit_cell)

            except Exception as e:
                print(f"Reflection path is: {reflections.path}")
                print(e)

        self.session.commit()

        print("Populated systems")
        print(
            (
                f"Got {self.session.query(func.count(Resolution.id)).scalar()} resolutions\n"
                f"Got {self.session.query(func.count(Spacegroup.id)).scalar()} spacegroups\n"
            )
        )

    def populate_panddas_errors(self, pandda_dirs_dir: Path):

        def tail(file: Path, n: int = 20) -> str:
            command: str = f"tail -n {n} {str(file)}"
            print(command)
            p = subprocess.Popen(command,
                                 shell=True,
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 )
            stdout, stderr = p.communicate()
            return str(stdout)

        for system in self.session.query(System).all():
            pandda_dir = pandda_dirs_dir / system.system

            pandda_json_file = pandda_dir / xlib.data.Constants.PANDDA_LOG_FILE

            # with open(pandda_json_file, "r") as f:
            #     pandda_json = json.load(f)

            # Determine if ran by looking for event table
            event_table_file = pandda_dir / xlib.data.Constants.PANDDA_ANALYSES_DIR / xlib.data.Constants.PANDDA_ANALYSE_EVENTS_FILE
            if event_table_file.exists():
                success = True
            else:
                success = False

            pandda = PanDDA(
                path=str(pandda_dir),
                system=system,
                success=success,
                runtime=0,
            )
            self.session.add(pandda)

            if not success:
                error_file = pandda_dirs_dir / xlib.data.Constants.PANDDA_JOB_ERROR_FILE.format(
                    system_name=system.system)
                pandda_error = PanDDAError(
                    path=str(error_file),
                    error=str(tail(Path(error_file))),
                    pandda=pandda,
                    system=system,
                )
                self.session.add(pandda_error)

        self.session.commit()

    def populate_autobuilds(self, autobuild_dirs_dir: Path):

        build_dict: xlib.data.BuildDict = xlib.data.BuildDict.from_autobuild_dir(autobuild_dirs_dir)

        for build_id in build_dict:
            build: xlib.data.Build = build_dict[build_id]

            system = self.session.query(System).filter(System.system == build_id.system.system).first()
            pandda = self.session.query(PanDDA).filter(PanDDA.system_id == system.id).first()
            dataset = self.session.query(Dataset).filter(Dataset.dtag == build_id.dtag.dtag).first()
            event = self.session.query(Event).filter(Event.system_id == system.id,
                                                     Event.dataset_id == dataset.id,
                                                     Event.event_idx == build_id.event_idx.event_idx
                                                     ).first()

            autobuild = Autobuild(path=str(build.build_file),
                                  system=system,
                                  dataset=dataset,
                                  pandda=pandda,
                                  event=event,
                                  )

            self.session.add(autobuild)

    def populate_events(self):

        for pandda in self.session.query(PanDDA):
            pandda_dir = Path(pandda.path)
            pandda_analyses_dir = pandda_dir / xlib.data.Constants.PANDDA_ANALYSES_DIR
            pandda_event_table_file = pandda_analyses_dir / xlib.data.Constants.PANDDA_ANALYSE_EVENTS_FILE

            # Check if processed
            if not pandda_event_table_file.exists():
                continue

            # Get table
            pandda_event_table: pd.DataFrame = pd.read_csv(str(pandda_event_table_file))
            #

            # Get events
            for index, row in pandda_event_table.iterrows():
                system = pandda.system
                dataset = self.session.query(Dataset).filter(Dataset.dtag == row["dtag"]).first()

                event_dir = pandda_dir / xlib.data.Constants.PANDDA_PROCESSED_DATASETS_DIR / dataset.dtag

                event_map = event_dir / xlib.data.Constants.PANDDA_EVENT_MAP_FILE.format(dtag=dataset.dtag,
                                                                                         event_idx=row["event_idx"],
                                                                                         bdc=round(row["bdc"], 2),
                                                                                         )

                event = Event(
                    event_idx=row["event_idx"],
                    bdc=row["1-BDC"],
                    x=row["x"],
                    y=row["y"],
                    z=row["z"],
                    analysed_Resolution=row["analysed_resolution"],
                    event_map=str(event_map),
                    system=system,
                    dataset=dataset,
                    pandda=pandda,
                )

                self.session.add(event)

        self.session.commit()

    def populate_test_train(
            self,
            database_path: str,
            out_dir: str,
            split: float,
    ):

        database_path = Path(database_path)

        # Get the systems
        system_list = self.session.query(System).all()

        # Get number of systems to lample
        num_systems = np.floor(len(system_list) * split)

        # Sample systems
        train_system_list = np.choice(np.arrange(0, len(system_list)), num_systems, replace=False)
        test_system_list = [x for x in np.arrange(0, len(system_list)) if x not in train_system_list]

        # Get the number of datasets for each
        train_dataset_list = self.session.query(Dataset).filter(Dataset.system_id.in_(train_system_list)).all()
        test_dataset_list = self.session.query(Dataset).filter(Dataset.system_id.in_(test_system_list)).all()
        print(
            (
                f"Train dataset: {len(train_dataset_list)}; test datasets: {len(test_dataset_list)}"
            )
        )

        # Get the number of reference models for each
        train_reference_list = self.session.query(ReferenceModel).filter(
            ReferenceModel.system_id.in_(train_system_list)).all()
        test_reference_list = self.session.query(ReferenceModel).filter(
            ReferenceModel.system_id.in_(test_system_list)).all()
        print(
            (
                f"Train reference models: {len(train_reference_list)}; test reference models: {len(test_reference_list)}"
            )
        )

        train_system_list = self.session.query(System).filter(System.id.in_(train_system_list)).all()
        test_system_list = self.session.query(System).filter(System.id.in_(test_system_list)).all()

        for system in train_system_list:
            train = Train(system=system)
            self.session.add(train)

        for system in test_system_list:
            test = Test(system=system)
            self.session.add(test)

        self.session.commit()

    # def populate_annotations(self, cutoff=5.0):

    #     reference_model_query = self.session.query(ReferenceModel)

    #     for event in self.session.query(Event).all():
    #         reference_model = reference_model_query.filter(ReferenceModel.dataset_id == event.dataset_id).first()

    #         if not reference_model:
    #             ann= False
    #         else:
    #             reference_structure = 
    #             reference_ligands = 
    #             reference_ligand_distances = 
    #             min_reference_ligand_distance = 
    #             if min_reference_ligand_distance < cutoff:
    #                 ann = True
    #             else:
    #                 ann = False

    #         annotation = Annotation(event=,
    #                                 dataset=,)

    #         self.session.add(annotation)

    #     self.session.commit()


@dataclasses.dataclass()
class Args:
    system_dirs_dir: Path
    pandda_dirs_dir: Path
    autobuild_dirs_dir: Path
    database_file: Path
    reference_model_dir: Path
    sequence_identity_dir: Path
    split: float

    @staticmethod
    def from_cmd():
        fields: Tuple[dataclasses.Field, ...] = dataclasses.fields(Args)
        parser = argparse.ArgumentParser()
        for field in fields:
            parser.add_argument(f"--{field.name}")
        args = parser.parse_args()

        args_dict = vars(args)

        typed_args = {field.name: field.type(args_dict[field.name]) for field in fields}

        return Args(**typed_args)


def main():
    args = Args.from_cmd()
    database = Database(
        args.database_file,
        overwrite=True,
    )

    database.populate_systems(args.system_dirs_dir)
    database.populate_sequence_identity(args.sequence_identity_dir)

    database.populate_models_reflections_compounds_smiles_datasets()
    database.populate_reference_models(args.reference_model_dir)

    database.populate_panddas_errors(args.pandda_dirs_dir)
    database.populate_events()

    database.populate_autobuilds(args.autobuild_dirs_dir)

    database.populate_resolution_spacegroup_unit_cell()  # long

    database.populate_test_train(args.split)

    # database.populate_annotations()

    # database.populate_autobuild_rmsds()
    # database.populate_autobuild_rsccs()
    # database.populate_autobuild_skeleton_scores()
    # database.populate_human_autobuild_comparison()
    # database.populate_real_space_clusterings()
    # database.populate_event_scores()


def make_database(database_file: str):
    database_file = Path(database_file)
    database = Database(database_file, overwrite=True)


def update_systems(database_file: str, system_dirs_dir: str):
    database_file = Path(database_file)
    system_dirs_dir = Path(system_dirs_dir)
    database = Database(database_file)

    database.populate_systems(system_dirs_dir)


def update_datasets(database_file: str):
    database_file = Path(database_file)
    database = Database(database_file)

    database.populate_models_reflections_compounds_smiles_datasets()
    database.populate_resolution_spacegroup_unit_cell()  # long


def update_panddas(database_file: str, pandda_dirs_dir: str):
    database_file = Path(database_file)
    pandda_dirs_dir = Path(pandda_dirs_dir)
    database = Database(database_file)

    database.populate_panddas_errors(pandda_dirs_dir)
    database.populate_events()


def update_clusters(database_file: str, cluster_dirs_dir: str):
    database_file = Path(database_file)
    cluster_dirs_dir = Path(cluster_dirs_dir)
    database = Database(database_file)

    database.populate_clusters(cluster_dirs_dir)


def update_autobuilds(database_file: str, autobuild_dirs_dir: str):
    database_file = Path(database_file)
    autobuild_dirs_dir = Path(autobuild_dirs_dir)
    database = Database(database_file)

    database.populate_autobuilds(autobuild_dirs_dir)


if __name__ == "__main__":
    # If run as a script creates and populates database
    fire.Fire(
        {
            "make": make_database,
            "update_systems": update_systems,
            "update_datasets": update_datasets,
            "update_panddas": update_panddas,
            "update_clusters": update_clusters,
            "update_autobuilds": update_autobuilds,
            "update_all": main,
        }
    )
