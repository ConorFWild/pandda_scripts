import os
from pathlib import Path
import json

import argparse
import dataclasses
from typing import *

from sqlalchemy import Column, ForeignKey, Integer, String, Float, Boolean, create_engine
from sqlalchemy.orm import relationship, sessionmaker
from sqlalchemy.ext.declarative import declarative_base

import pandas as pd

import xlib


class Constants:
    COMPOUND_TABLE = "compound"
    REFLECTIONS_TABLE = "reflections"
    MODEL_TABLE = "model"
    DATASET_TABLE = "dataset"
    SYSTEM_TABLE = "system"
    PANDDA_TABLE = "pandda"
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
    
    

base = declarative_base()
        
class Compound(base):
    __tablename__ = Constants.COMPOUND_TABLE
    id = Column(Integer, primary_key=True)
    path = Column(String(255))

class Reflections(base):
    __tablename__ = Constants.REFLECTIONS_TABLE
    id = Column(Integer, primary_key=True)
    path = Column(String(255))

class Model(base):
    __tablename__ = Constants.MODEL_TABLE
    id = Column(Integer, primary_key=True)
    path = Column(String(255))

class Dataset(base):
    __tablename__ = Constants.DATASET_TABLE
    id = Column(Integer, primary_key=True)
    dtag = Column(String(255))
    path = Column(String(255))
    
    # Foriegn keys
    reflections_id = Column(Integer, ForeignKey(Reflections.id))
    model_id = Column(Integer, ForeignKey(Model.id))
    compound_id = Column(Integer, ForeignKey(Compound.id))
    
    # Relationships
    reflections = relationship(Reflections)
    model = relationship(Model)
    compound = relationship(Compound)
    

class System(base):
    __tablename__ = Constants.SYSTEM_TABLE
    id = Column(Integer, primary_key=True)
    system = Column(String(255))
    path = Column(String(255))

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
    

class Event(base):
    __tablename__ = Constants.EVENT_TABLE
    id = Column(Integer, primary_key=True)
    event_idx = Column(Integer)
    x = Column(Float)
    y = Column(Float)
    z= Column(Float)
    analysed_Resolution = Column(Float)
    bdc = Column(Float)
    
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
    dataset_id= Column(Integer, ForeignKey(Dataset.id))
    
    # Relationships
    system = relationship(System)
    dataset = relationship(Dataset)
    
    
class Autobuild(base):
    __tablename__ = Constants.AUTOBUILD_TABLE
    id = Column(Integer, primary_key=True)
    path = Column(String(255))
    
    # Foreign keys
    system_id = Column(Integer, ForeignKey(System.id))
    dataset_id= Column(Integer, ForeignKey(Dataset.id))
    pandda_id= Column(Integer, ForeignKey(PanDDA.id))
    event_id= Column(Integer, ForeignKey(Event.id))
    
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
    dataset_id= Column(Integer, ForeignKey(Dataset.id))
    pandda_id= Column(Integer, ForeignKey(PanDDA.id))
    event_id= Column(Integer, ForeignKey(Event.id))
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
    dataset_id= Column(Integer, ForeignKey(Dataset.id))
    pandda_id= Column(Integer, ForeignKey(PanDDA.id))
    event_id= Column(Integer, ForeignKey(Event.id))
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
    dataset_id= Column(Integer, ForeignKey(Dataset.id))
    pandda_id= Column(Integer, ForeignKey(PanDDA.id))
    event_id= Column(Integer, ForeignKey(Event.id))
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
    dataset_id= Column(Integer, ForeignKey(Dataset.id))
    pandda_id= Column(Integer, ForeignKey(PanDDA.id))
    event_id= Column(Integer, ForeignKey(Event.id))
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
    dataset_id= Column(Integer, ForeignKey(Dataset.id))
    pandda_id= Column(Integer, ForeignKey(PanDDA.id))
    event_id= Column(Integer, ForeignKey(Event.id))
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
    dataset_id= Column(Integer, ForeignKey(Dataset.id))
    pandda_id= Column(Integer, ForeignKey(PanDDA.id))
    event_id= Column(Integer, ForeignKey(Event.id))

    # Relationships
    system = relationship(System)
    dataset = relationship(Dataset)
    pandda = relationship(PanDDA)
    event = relationship(Event)

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
            system_path: str =str(system_path) 
        
            system = System(
                system=system_name, 
                path=system_path,
                )
            self.session.add(system)
        
        self.session.commit()
        
    def populate_models_reflections_compounds_datasets(self):
        
        for system in self.session.query(System):
            system_path = Path(system.path)
            
            dataset_dirs_dir = system_path / system.system
            
            dataset_dir_list = list(path for path in dataset_dirs_dir.glob("*") if path.is_dir())
            
            for dataset_dir in dataset_dir_list:
                compound_path = dataset_dir / Constants.COMPOUND_FILE
                reflections_path = dataset_dir / Constants.REFLECTIONS_FILE
                model_path = dataset_dir / Constants.MODEL_FILE
            
                model = Model(path=str(compound_path))
            
                reflections = Reflections(path=str(reflections_path))
                
                compound = Compound(path=str(model_path))
                
                dataset = Dataset(dtag=dataset_dir.name,
                                  path=str(dataset_dir),
                                  reflections=reflections,
                                  compound=compound,
                                  model=model,
                                  )
                
                self.session.add(model)
                self.session.add(reflections)
                self.session.add(compound)
                self.session.add(dataset)
                
        self.session.commit()
        
    def populate_panddas(self, pandda_dirs_dir: Path):
        
        for system in self.session.query(System):
            pandda_dir = pandda_dirs_dir / system.system
            
            pandda_json_file = pandda_dir / xlib.Constants.PANDDA_LOG_FILE
            
            with open(pandda_json_file, "r") as f:
                pandda_json = json.load(f)
        
            pandda = PanDDA(
                path=str(pandda_dir),
                system=system,
                success=pandda_json["success"],
                runtime=pandda_json["runtime"],
                )
            self.session.add(pandda)
        
        self.session.commit()
        
    def populate_autobuilds(self, autobuild_dirs_dir: Path):
        
        build_dict: xlib.BuildDict = xlib.BuildDict.from_autobuild_dir(autobuild_dirs_dir)
        
        for build_id in build_dict:
            build: xlib.Build = build_dict[build_id]
            
            system = self.session.query(System).filter(System.system == build_id.system.system).first()
            pandda = self.session.query(PanDDA).filter(PanDDA.system == build_id.system.system).first()
            dataset = self.session.query(Dataset).filter(Dataset.dtag == build_id.dtag.dtag ).first()
            event = self.session.query(Event).filter(Event.system.system == build_id.system.system &
                                                     Event.dataset.dtag == build_id.dtag.dtag &
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
            pandda_analyses_dir = pandda_dir / xlib.Constants.PANDDA_ANALYSES_DIR
            pandda_event_table_file = pandda_analyses_dir / xlib.Constants.PANDDA_ANALYSE_EVENTS_FILE
            
            # Get table
            pandda_event_table: pd.DataFrame = pd.read_csv(str(pandda_event_table_file))
            
            # Get events
            for index, row in pandda_event_table.iterrows():
                system = pandda.system
                dataset = self.session.query(Dataset).filter(Dataset.dtag == row["dtag"]).first()
                
                event_dir = pandda_dir/ xlib.Constants.PANDDA_PROCESSED_DATASETS_DIR / dataset.dtag
                
                
                event = Event(
                    event_idx = row["event_idx"],
                    bdc = row["1-BDC"],
                    x = row["x"],
                    y = row["y"],
                    z = row["z"],
                    analysed_Resolution = row["analysed_resolution"],
                    system=system,
                    dataset=dataset,
                    pandda=pandda,
                    )
                            
                self.session.add(event)
        
        self.session.commit()
        
        
@dataclasses.dataclass()
class Args:
    system_dirs_dir: Path
    pandda_dirs_dir: Path
    autobuild_dirs_dir: Path
    database_file: Path
    
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
    database.populate_models_reflections_compounds_datasets()
    database.populate_panddas(args.pandda_dirs_dir)
    database.populate_autobuilds(args.autobuild_dirs_dir)
    database.populate_events()
    
    database.populate_reference_models()
    database.populate_autobuild_rmsds()
    database.populate_autobuild_rsccs()
    database.populate_autobuild_skeleton_scores()
    database.populate_human_autobuild_comparison()
    database.populate_real_space_clusterings()
    database.populate_event_scores()
    
if __name__ == "__main__":
    # If run as a script creates and populates database
    main()