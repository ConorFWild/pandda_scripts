from __future__ import annotations

import os

from typing import *
import dataclasses

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

import tables

import xlib

# ############
# Cosntants
# ############
class TableConstants:
    # Groups
    SYSTEM_GROUP_NAME = "system"
    PANDDA_GROUP_NAME = "pandda"
    EVENT_GROUP_NAME = "event"
    AUTOBUILD_GROUP_NAME = "autobuild"
    
    # Tables
    SYSTEM_TABLE_NAME = "system"
    PANDDA_TABLE_NAME = "pandda"
    EVENT_TABLE_NAME = "event"
    AUTOBUILD_TABLE_NAME = "autobuild"
    
    # Defaults
    DEFAULT_DATABASE_FILE = "/data/share-2/conor/pandda/output/database.hdf5"
    DEFAULT_SYSTEM_DIRS_DIR = "/data/share-2/conor/pandda/data/pandda_inputs"
    DEFAULT_PANDDA_DIRS_DIR = "/data/share-2/conor/pandda/output/all_systems"
    DEFAULT_AUTOBUILD_DIRS_DIR = "/data/share-2/conor/pandda/output/all_autobuild"
    
# #############
# Records
# ##############

class SystemRecord(tables.IsDescription):
    system = tables.StringCol(255)
    path = tables.StringCol(255)

@dataclasses.dataclass()
class BuildRecord():
    @staticmethod
    def from_row(row) -> PanDDARecord:
        return EventRecord(
            row["system"],
            row["dtag"] ,
            row["event_idx"],
            row["build_cluster"] ,
            row["build_number"] ,
            row["rscc"] ,
            row["file"],
        )

    def fill_row(self, row: tables.tableextension.Row) -> tables.tableextension.Row:
        row["system"] = self.system
        row["dtag"] = self.dtag
        row["event_idx"] = self.event_idx
        row["build_cluster"] = self.build_cluster_id
        row["build_number"] = self.build_number_id
        row["rscc"] = self.build_rscc
        row["file"] = self.build_file
        return row

class BuildRecordDescription(tables.IsDescription):
    system = tables.StringCol(255)
    dtag = tables.StringCol(255)
    event_idx = tables.Int32Col()
    build_cluster = tables.Int32Col()
    build_number = tables.Int32Col()
    rscc = tables.Float32Col()
    file = tables.StringCol(255)

@dataclasses.dataclass()
class PanDDARecord:
    system: str
    path: str
    event_table_file: str

    @staticmethod
    def from_row(row) -> PanDDARecord:
        return PanDDARecord(
            row["system"],
            row["path"] ,
            row["event_table_file"],
        )

    def fill_row(self, row: tables.tableextension.Row) -> tables.tableextension.Row:
        row["system"] = self.system
        row["path"] = self.pandda_dir
        row["event_table_file"] = self.event_table_file
        return row
        
    
class PanDDARecordDescription(tables.IsDescription):
    system = tables.StringCol(255)
    path = tables.StringCol(255)
    event_table_file = tables.StringCol(255)

class MTZRecord(tables.IsDescription):
    system = tables.StringCol(255)
    dtag = tables.StringCol(255)
    path = tables.StringCol(255)

@dataclasses.dataclass()
class EventRecord:
    system: str
    path: str
    event_idx: int
    x: float
    y: float
    z: float
    bdc: float
    resolution: float
    
    @staticmethod
    def from_pandda_event_table_row(row, system) -> PanDDARecord:
        return EventRecord(
            system,
            row["dtag"] ,
            row["event_idx"],
            row["x"] ,
            row["y"] ,
            row["z"] ,
            row["1-BDC"],
            row["analysed_resolution"],
        )        
    
    @staticmethod
    def from_row(row) -> PanDDARecord:
        return EventRecord(
            row["system"],
            row["dtag"] ,
            row["event_idx"],
            row["x"] ,
            row["y"] ,
            row["z"] ,
            row["bdc"],
            row["resolution"],
        )

    def fill_row(self, row: tables.tableextension.Row) -> tables.tableextension.Row:
        row["system"] = self.system
        row["dtag"] = self.dtag
        row["event_idx"] = self.event_idx
        row["x"] = self.x
        row["y"] = self.y
        row["z"] = self.z
        row["bdc"] = self.bdc
        row["resolution"] = self.resolution
        return row

class EventRecordDescription(tables.IsDescription):
    system = tables.StringCol(255)
    dtag = tables.StringCol(255)
    event_idx = tables.Int32Col()
    x = tables.Float32Col()
    y = tables.Float32Col()
    z = tables.Float32Col()
    bdc = tables.Float32Col()
    resolution = tables.Float32Col()

class StructureRecord(tables.IsDescription):
    system = tables.StringCol(255)
    dtag = tables.StringCol(255)
    path = tables.StringCol(255)

@dataclasses.dataclass()
class Database:
    _file: Path 
    _table: tables.File
    
    system_group: tables.Group
    pandda_group: tables.Group
    event_group: tables.Group
    autobuild_group: tables.Group
    
    system_table: tables.Table
    pandda_table: tables.Table
    event_table: tables.Table
    autobuild_table: tables.Table
                        
    
    @staticmethod
    def from_file(path: str = TableConstants.DEFAULT_DATABASE_FILE, overwrite: bool=True) -> Database:
        
        _file: Path = Path(path)
        
        if overwrite:
            os.remove(str(_file))
        
        _table = tables.open_file(
            str(path), 
            "a",
            )
        
        # Create groups
        system_group: tables.Group = _table.create_group(_table.root, TableConstants.SYSTEM_GROUP_NAME)
        pandda_group: tables.Group = _table.create_group(_table.root, TableConstants.PANDDA_GROUP_NAME)
        event_group: tables.Group = _table.create_group(_table.root, TableConstants.EVENT_GROUP_NAME)
        autobuild_group: tables.Group = _table.create_group(_table.root, TableConstants.AUTOBUILD_GROUP_NAME)
        
        # Create tables
        system_table = _table.create_table(system_group, 
                            TableConstants.SYSTEM_TABLE_NAME,
                            SystemRecord,
                            )
        pandda_table = _table.create_table(pandda_group, 
                            TableConstants.PANDDA_TABLE_NAME,
                            PanDDARecordDescription,
                            )
        event_table = _table.create_table(event_group, 
                            TableConstants.EVENT_TABLE_NAME,
                            EventRecordDescription,
                            )
        autobuild_table = _table.create_table(autobuild_group, 
                            TableConstants.AUTOBUILD_TABLE_NAME,
                            BuildRecordDescription,
                            )
        
        return Database(_file, _table,
                        system_group,
                        pandda_group,
                        event_group,
                        autobuild_group,
                        system_table,
                        pandda_table,
                        event_table,
                        autobuild_table,               
                        )
        
    def make(self) -> None:
        # PanDDAs
        self._table.create_group(
            self._table.root,
            TableConstants.PANDDA_RESULTS,
            )
        
        self._table.create_table(
            self._table.root,
            
            )
        
        # Autobuilds
        
    def populate_systems(self, path: Path) -> None:
        system_table: tables.Table = self.system_table
        row: tables.tableextension.Row = system_table.row
        
        system_path_list: List[Path] = list(path for path in path.glob("*") if path.is_dir())
        
        for system_path in system_path_list:
            system: str = system_path.name
            system_path: str =str(system_path) 
            
            # Create record
            row = Database.fill_row_system(
                row,
                system,
                system_path,
                )
            row.append()
            
        self._table.flush()
        print(system_table)
        
    def populate_panddas(self, pandda_dirs_dir: Path) -> None:
        system_table: tables.Table = self.system_table
        pandda_table: tables.Table = self.pandda_table
        pandda_row: tables.tableextension.Row = pandda_table.row
        
        # Get event tables
        for system_record in system_table.iterrows():
            system: str = str(system_record["system"], 'utf-8')
            
            event_table_file: Path = pandda_dirs_dir / system / xlib.Constants.PANDDA_ANALYSES_DIR / xlib.Constants.PANDDA_ANALYSE_EVENTS_FILE
            print(event_table_file)
            if event_table_file.exists():
                pandda_dir: Path = pandda_dirs_dir / system
                    
                # Make record
                pandda_row = Database.fill_row_pandda(
                    pandda_row,
                    system,
                    pandda_dir,
                    event_table_file,
                )
                pandda_row.append()
                
        # Flush records
        self._table.flush()
        print(pandda_table)
        
    def populate_events(self, pandda_dirs_dir, autobuild_dirs_dir) -> None:
        system_table: tables.Table = self.system_table
        pandda_table: tables.Table = self.pandda_table
        event_table: tables.Table = self.event_table
        event_row:  tables.tableextension.Row = event_table.row
        
        for pandda_row in pandda_table:
            pandda_record: PanDDARecord = PanDDARecord.from_row(pandda_row)
            event_table_file = str(pandda_record.event_table_file, "utf-8")
            # Get table
            pandda_event_table: pd.DataFrame = pd.read_csv(str(event_table_file))
            
            # Get events
            for index, row in pandda_event_table.iterrows():
                event_record: EventRecord = EventRecord.from_pandda_event_table_row(row,
                                                                                    pandda_record.system,
                                                                                    )
                
                event_record.fill_row(event_row)
                
                event_row.append()
                
            # Flush records
            self._table.flush()
            print(event_table)
            

        
    def populate_autobuilds(self, autobuild_dirs_dir: Path) -> None:
        build_table: tables.Table = self.event_table
        build_row:  tables.tableextension.Row = build_table.row
        
        build_dict: xlib.BuildDict = xlib.BuildDict.from_autobuild_dir(autobuild_dirs_dir)
        
        for build_id in build_dict:
            build: xlib.Build = build_dict[build_id]

            # Get event record
            build_row.system = build_id.system.system
            build_row.dtag = build_id.dtag.dtag
            build_row.event_idx = build_id.event_idx.event_idx
            build_row.build_cluster = build_id.build_cluster.build_cluster_id
            build_row.build_number = build_id.build_number.build_number_id
            build_row.rscc = build.build_rscc
            build_row.file = build.build_file

        
        self._table.flush()
        print(build_table)
        
    @staticmethod
    def fill_row_system(
        row: tables.tableextension.Row,
        system: str,
        system_path: str,
        ) -> tables.tableextension.Row:
        row["system"] = system
        row["path"] = system_path
        return row
    
    @staticmethod
    def fill_row_pandda(
        pandda_row: tables.tableextension.Row,
        system: str,
        pandda_dir: str,
        event_table_file: str,
        
    ) -> tables.tableextension.Row:
        pandda_row["system"] = system
        pandda_row["path"] = pandda_dir
        pandda_row["event_table_file"] = event_table_file
        return pandda_row

# ############
# Args
# ##############
@dataclasses.dataclass()
class Args:
    database_file: Path
    system_dirs_dir: Path
    pandda_dirs_dir: Path
    autobuild_dirs_dir: Path
    debug: int
    
    @staticmethod
    def from_args(args: Any):
        database_file = Path(args.database_file)
        system_dirs_dir = Path(args.system_dirs_dir)
        pandda_dirs_dir = Path(args.pandda_dirs_dir)
        autobuild_dirs_dir = Path(args.autobuild_dirs_dir)
        debug: int = int(args.debug)
        
        return Args(
            database_file,
            system_dirs_dir,
            pandda_dirs_dir,
            autobuild_dirs_dir,
            debug,
        )
    
    @staticmethod
    def from_cmd():
        parser = argparse.ArgumentParser()
        parser.add_argument("--database_file",
                            default=TableConstants.DEFAULT_DATABASE_FILE,
                            )
        parser.add_argument("--system_dirs_dir",
                            default=TableConstants.DEFAULT_SYSTEM_DIRS_DIR,
                            )
        parser.add_argument("--pandda_dirs_dir",
                            default=TableConstants.DEFAULT_PANDDA_DIRS_DIR,
                            )
        parser.add_argument("--autobuild_dirs_dir",
                            default=TableConstants.DEFAULT_AUTOBUILD_DIRS_DIR,
                            )
        parser.add_argument("--debug",
                            default=2,
                            )
        args = parser.parse_args()
        
        return Args.from_args(args)


# ################
# Main
# #############
def main():
    # Get args
    args: Args = Args.from_cmd()
    
    # Create database - use default location
    database: Database = Database.from_file(args.database_file)

    # Populate systems
    database.populate_systems(args.system_dirs_dir)
    
    # populate panddas
    database.populate_panddas(args.pandda_dirs_dir)
    
    # Populate events
    database.populate_events(args.pandda_dirs_dir,
                             args.autobuild_dirs_dir,
                             )
    
    # Populate autobuilds
    database.populate_autobuilds(args.autobuild_dirs_dir)
    
    
    
# ##########
# Script
# ############
if __name__ == "__main__":
    main()