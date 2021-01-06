from __future__ import annotations

import os
from typing import Dict
import time
import psutil
import pickle
from shlex import split
from pprint import PrettyPrinter
from pathlib import Path
import json

import numpy as np
import pandas as pd

from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from hdbscan import HDBSCAN


import gemmi


from xlib.pandda_types import (JoblibMapper, PanDDAFSModel, Datasets, Reference, 
                                    Grid, Alignments, Shells, Xmaps, Xmap,
                                    XmapArray, Model, Dtag, Zmaps, Clusterings,
                                    Events, SiteTable, EventTable,
                                    JoblibMapper, Event
                                    )
import joblib
from joblib.externals.loky import set_loky_pickler
set_loky_pickler('pickle')

import plotly.graph_objects as go

        
@dataclasses.dataclass()
class Args:
    system_dirs_dir: Path
    out_dir: Path
    
    @staticmethod
    def from_cmd():
        
        fields: Tuple[dataclasses.Field, ...] = dataclasses.fields(Args)
        parser = argparse.ArgumentParser()
        for field in fields:
            parser.add_argument(f"--{field.name}")
        args = parser.parse_args()
        args_dict = vars(args)
        typed_args = [field.type(args_dict[field.name]) for field in fields]
        
        return Args.from_args(args)
    
    
class Constants:
    RESIDUE_CLUSTER_PLOT_FILE = ""
    CLUSTERINGS_PLOT_FILE = ""
    RECORDS_JSON_FILE = ""


def select_partition(xmap: Xmap, partition):
    array = xmap.to_array()
    selection = (
        np.array(partition_coord[0] for partition_coord in partition),
        np.array(partition_coord[1] for partition_coord in partition),
        np.array(partition_coord[2] for partition_coord in partition)
                 )
    
    selected_values_array = array[selection]
    
    return selected_values_array

def make_sample_by_feature_matrix(selection_dict):
    sample_by_feature_array = np.vstack([selection for selection in selection_dict.values()])
    
    return sample_by_feature_array

        
def embed_pca(sample_by_feature_matrix):
    pca_obj = PCA(n_components=50)
    
    pca_obj.fit_transform(sample_by_feature_matrix)


def embed_tsne(pca_matrix):
    tsne_obj = TSNE(n_components=2)
    
    tsne_obj.fit_transform(pca_matrix)

def cluster(tsne_matrix):
    hdbscan_obj = HDBSCAN()
    
    hdbscan_obj.fit(tsne_matrix)
    
    return hdbscan_obj.labels_


def make_clustering_record(selection_dict, tsne_matrix, clustering): 
    records = {}
    for dtag, tsne_coords, clustering in zip(selection_dict, tsne_matrix, clustering):
        records[dtag] = {
            "x": tsne_coords[0],
            "y": tsne_coords[1],
            "cluster": clustering,
        }
        
    return records

def save_records(records, out_dir):
    with open(out_dir / Constants.RECORDS_JSON_FILE, "w") as f:
        json.dump(records, f)

def plot_clustering(residue_id, tsne_matrix, clustering, out_dir):
    fig = go.Figure()

    
    fig.add_trace(
        go.Scatter(
            x=tsne_matrix[:,0],
            y=tsne_matrix[:,1],
            mode='markers',
            name='Original Plot',
            marker=dict(
                color=clustering,
            ),
        )
    )
    
    fig.write_html(str(out_dir / Constants.CLUSTERINGS_PLOT_FILE))


def plot_clusterings(records, out_dir: Path):

    fig = go.Figure()
    
    # Invert dict
    clustering_dict_dict = {}
    for residue_id in records:
        residue_record = records[residue_id]
        for dtag in residue_record:
            if dtag.dtag not in clustering_dict_dict:
                clustering_dict_dict[dtag.dtag] = {}
            
            clustering_dict_dict[dtag.dtag][residue_id] = residue_record[dtag]

    # Add traces to graph
    for dtag in clustering_dict_dict:
        series = pd.Series.from_dict(clustering_dict_dict[dtag])
        
        fig.add_trace(go.Scatter(
            y=series,
            mode='lines+markers',
            name='Original Plot'
        ))
    
    fig.write_html(str(out_dir / Constants.CLUSTERINGS_PLOT_FILE))


def main():
    ###################################################################
    # # Configuration
    ###################################################################
    print("Getting config")
    args = Args.from_cmd()
    
    print("Getting multiprocessor")
    mapper = JoblibMapper.initialise()
    
    
    ###################################################################
    # # Pre-pandda
    ###################################################################
    
    # Get datasets
    print("Loading datasets")
    datasets_initial: Datasets = Datasets.from_dir(args.data_dirs)
    
    # datasets_initial: Datasets = datasets_initial.trunate_num_datasets(100)

    # Initial filters
    print("Filtering invalid datasaets")
    datasets_invalid: Datasets = datasets_initial.remove_invalid_structure_factor_datasets(args.structure_factors)

    datasets_low_res: Datasets = datasets_invalid.remove_low_resolution_datasets(
        args.low_resolution_completeness)
    
    datasets_rfree: Datasets = datasets_low_res.remove_bad_rfree(args.filtering.max_rfree)
    
    datasets_wilson: Datasets = datasets_rfree.remove_bad_wilson(args.max_wilson_plot_z_score)  # TODO
    
    # Select refernce
    print("Getting reference")
    reference: Reference = Reference.from_datasets(datasets_wilson)
    
    # Post-reference filters
    print("smoothing")
    start = time.time()
    datasets_smoother: Datasets = datasets_wilson.smooth_datasets(reference, 
                                                structure_factors=args.structure_factors,
                                                mapper=mapper,
                                                )
    finish = time.time()
    print(f"Smoothed in {finish-start}")  

    print("Removing dissimilar models")
    datasets_diss_struc: Datasets = datasets_smoother.remove_dissimilar_models(reference,
                                                        args.max_rmsd_to_reference,
                                                        )
    
    datasets_diss_space: Datasets = datasets_diss_struc.remove_dissimilar_space_groups(reference)

    datasets = datasets_diss_space

    # Grid
    print("Getting grid")
    grid: Grid = Grid.from_reference(reference,
                            args.outer_mask,
                                args.inner_mask_symmetry,
                                    sample_rate=args.sample_rate,
                                )

    print("Getting alignments")
    alignments: Alignments = Alignments.from_datasets(
        reference,
        datasets,
        )
    
        
    ###################################################################
    # # Process shells
    ###################################################################
    sorted_dtags = list(sorted(datasets.datasets.keys(),
                            key=lambda dtag: datasets[dtag].reflections.resolution().resolution,
                            ))
    
    min_res = datasets[sorted_dtags[-1]].reflections.resolution().resolution

    
    print("Truncating datasets")
    truncated_datasets: Datasets = datasets.truncate(resolution=min_res,
                                                                structure_factors=args.structure_factors,
                                                                )

    # Get maps
    xmaps = Xmaps.from_aligned_datasets_c(
        truncated_datasets, 
        alignments, 
        grid,
        args.structure_factors, 
        sample_rate=args.sample_rate,
        mapper=mapper,
        ) # n x (grid size) with total_mask > 0
    
    records = {}
    # Iterate over residues
    for residue_id in reference.dataset.structure.protein_residue_ids():
        partition = grid.partitioning[residue_id]
        
        selection_dict = {dtag: select_partition(xmap, partition) for dtag, xmap in xmaps.xmaps.items()}
        
        sample_by_feature_matrix = make_sample_by_feature_matrix(selection_dict)
        
        pca_matrix = embed_pca(sample_by_feature_matrix)
        
        tsne_matrix = embed_tsne(pca_matrix)
        
        clustering = cluster(tsne_matrix)
        
        plot_clustering(tsne_matrix, clustering, args.out_dir)
        
        record = make_clustering_record(selection_dict, tsne_matrix, clustering)
        
        records[residue_id] = record


    plot_clusterings(records)
    
    save_records(records, args.out_dir)
    
    
    
        
        
