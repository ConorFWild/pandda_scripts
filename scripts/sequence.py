import os
from pathlib import Path
import json

import argparse
import dataclasses
from typing import *
import subprocess

import pandas as pd

import gemmi

import pandas as pd

import xlib

from xlib import data
from xlib import database_sql


@dataclasses.dataclass()
class Args:
    database_file: Path
    identity_table_dir: Path
    
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
    database = database_sql.Database(
        args.database_file, 
        overwrite=False,
        )
    
    
    # Get models
    example_model_dict = {}
    for system in database.session.query(database_sql.System).all():
        print(f"Getting model for system: {system.system}")
        example_dataset = database.session.query(
            database_sql.Dataset
            ).filter(
                database_sql.Dataset.system_id == system.id,
                ).first()

        model = example_dataset.model
        
        structure_path = model.path
        
        structure = gemmi.read_structure(str(structure_path))
        
        structure.setup_entities()
        structure.assign_label_seq_id()
        
        example_model_dict[system.id] = structure
        
    print(f"Got {len(example_model_dict)} systems")
        
    # Get identity matrix
    record_list = []
    blosum62 = gemmi.prepare_blosum62_scoring()

    for reference_system, reference_structure in example_model_dict.items():
        record = {}
        record["system"] = reference_system
        for moving_system, moving_structure in example_model_dict.items():
            reference_sequence = gemmi.one_letter_code(reference_structure.entities[0].full_sequence)
            moving_sequence = gemmi.one_letter_code(moving_structure.entities[0].full_sequence)
            
            
            alignment = gemmi.align_string_sequences(
                list(reference_sequence),
                list(moving_sequence),
                [], blosum62)
            
            identity = alignment.calculate_identity()
            record[moving_system] = identity
            
        record_list.append(record)            
    
    print(len(record_list))

    table = pd.DataFrame(record_list)
    
    table.to_csv(str(args.identity_table_dir / data.Constants.IDENTITY_TABLE_FILE))
            
    print(table)
            
if __name__ == "__main__":
    main()