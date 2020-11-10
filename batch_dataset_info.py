from __future__ import annotations

from typing import Dict, List
from dataclasses import dataclass


from pathlib import Path

dataclass()
class Systems:
    systems: Dict[SystemID, System]
    
class System:
    datasets: Dict[DatasetID, Dataset]
    
class Dataset:
    structure: Structure
    reflections: Reflections
    ligand: Ligand

class Structure:
    path: Path
    rfree: RFree

class 

def analyse(path: Path):
    
    # Get number of datasets
    
    # Get resolution distribution
    
    # Get rfactor distribution
    

def main():
    # Get data dirs
    
    # Analyse each
    
    # Summarise

if __name__ == "__main__":
    main()