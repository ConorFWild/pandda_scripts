from __future__ import annotations

import os
from typing import *
import time

import psutil
import pickle
from shlex import split
from pprint import PrettyPrinter
from pathlib import Path
import json
import dataclasses
import argparse


import numpy as np
import pandas as pd

from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from hdbscan import HDBSCAN


import gemmi

from pandda_gemmi import pandda_types

from pandda_gemmi.pandda_types import (JoblibMapper, PanDDAFSModel, Datasets, Reference, Dataset,
                                    Grid, Alignments, Shells, Xmaps, Xmap,
                                    XmapArray, Model, Dtag, Zmaps, Clusterings,
                                    Events, SiteTable, EventTable,
                                    JoblibMapper, Event, ResidueID, StructureFactors,
                                       sample_residue
                                    )
import joblib
from joblib.externals.loky import set_loky_pickler
set_loky_pickler('pickle')

import plotly.graph_objects as go

from plotly import figure_factory as ff

from sklearn import decomposition
from sklearn import manifold

# from scipy import stats
from numpy import random
from sklearn import mixture
from sklearn.neighbors import DistanceMetric

import matplotlib
from matplotlib import pyplot as plt
from scipy.cluster import hierarchy 
import scipy.cluster.hierarchy as spc
from scipy.spatial import distance

from numpy import random
from sklearn import mixture
from sklearn.neighbors import DistanceMetric



def sample_residue(truncated_dataset: Dataset,
                   grid: Grid,       
                   residue_id,
                    alignment: Alignment, 
                    structure_factors: StructureFactors, 
                    sample_rate: float, 
                    ) -> List[float]:
        
    point_position_dict = grid.partitioning[residue_id]
    
    unaligned_xmap: gemmi.FloatGrid = truncated_dataset.reflections.reflections.transform_f_phi_to_map(structure_factors.f,
                                                                                                structure_factors.phi,
                                                                                                sample_rate=sample_rate,
                                                                                                )    
    
    unaligned_xmap_array = np.array(unaligned_xmap, copy=False)


    std = np.std(unaligned_xmap_array)
    unaligned_xmap_array[:, :, :] = unaligned_xmap_array[:, :, :] / std

    # Unpack the points, poitions and transforms
    point_list: List[Tuple[int, int, int]] = []
    position_list: List[Tuple[float, float, float]] = []
    transform_list: List[gemmi.transform] = []
    com_moving_list: List[np.array] = []
    com_reference_list: List[np.array] = []
            
    al = alignment[residue_id]
    transform = al.transform.inverse()
    com_moving = al.com_moving
    com_reference = al.com_reference
    
    for point, position in point_position_dict.items():
        
        point_list.append(point)
        position_list.append(position)
        transform_list.append(transform)
        com_moving_list.append(com_moving)
        com_reference_list.append(com_reference)
    
    sampled_points = gemmi.interpolate_to_list(unaligned_xmap,
                                               grid.grid,
                                    point_list,
                                 position_list,
                                 transform_list,
                                 com_moving_list,
                                 com_reference_list,           
                              )
    
    return np.array(sampled_points)


def embed(distance_matrix):
    pca = decomposition.PCA(n_components=50)
    tsne = manifold.TSNE(n_components=2)
    transform = pca.fit_transform(distance_matrix)
    transform = tsne.fit_transform(transform)
    return transform


@dataclasses.dataclass()
class Args:
    data_dirs: Path
    out_dir: Path
    pdb_regex: str
    mtz_regex: str
    structure_factors: StructureFactors = StructureFactors.from_string("2FOFCWT,PH2FOFCWT")
    # structure_factors: StructureFactors = StructureFactors.from_string("FWT,PHWT")
    low_resolution_completeness: float = 4.0
    max_rfree: float = 0.4
    max_wilson_plot_z_score: float = 5.0
    max_rmsd_to_reference: float = 1.5
    outer_mask: float = 6.0
    inner_mask_symmetry: float = 3.0
    sample_rate: float = 3.0
    num_components: int = 20
    
    @staticmethod
    def from_cmd():
        
        fields = dataclasses.fields(Args)
        parser = argparse.ArgumentParser()
        for field in fields:
            if field.default is None:
                parser.add_argument(f"--{field.name}")
        args = parser.parse_args()
        
        args_dict = vars(args)
        
        typed_args = [field.type(args_dict[field.name]) for field in fields if field.default is None]
        
        return Args(typed_args)
    
    @staticmethod
    def from_arg_list(arg_list):
        
        fields = dataclasses.fields(Args)
        parser = argparse.ArgumentParser()
        for field in fields:
            print(field.default)
            if type(field.default) is type(dataclasses.MISSING):
                print(field.name)
                parser.add_argument(f"--{field.name}")
        args = parser.parse_args(arg_list)
        
        args_dict = vars(args)
        
        print([field.type for field in fields])
        
        typed_args = [field.type(args_dict[field.name]) for field in fields if type(field.default) is type(dataclasses.MISSING)]
        
        return Args(typed_args)
    

def save_dendrogram_plot(linkage, 
                         labels, 
                         dendrogram_plot_file,
                         ):
    fig, ax = plt.subplots(figsize=(60,40))
    dn = spc.dendrogram(linkage, ax=ax, labels=labels, leaf_font_size=10)
    fig.savefig(str(dendrogram_plot_file))
    fig.clear()
    plt.close(fig)
    
    
def get_linkage_from_correlation_matrix(correlation_matrix):
    condensed = distance.squareform(1.0-correlation_matrix)
    linkage = spc.linkage(condensed, method='complete')

    return linkage


def save_correlation_plot(correlation_matrix, correlation_plot_file):
                    
    fig, ax = plt.subplots(figsize=(20,20))
    ax.imshow(correlation_matrix)
    fig.savefig(str(correlation_plot_file))
    fig.clear()
    plt.close(fig)


def cluster_linkage(linkage, cutoff):
    idx = spc.fcluster(linkage, cutoff, 'distance')

    return idx


def get_corr(reference_sample_mask, sample_mask, diag):
    reference_mask_size = np.sum(reference_sample_mask)
    sample_mask_size = np.sum(sample_mask)
    
    denominator = max(sample_mask_size, reference_mask_size)
    
    if denominator == 0.0:
        if diag:
            corr = 1.0
        else:
            corr = 0.0
    
    else: 

        corr = np.sum(sample_mask[reference_sample_mask == 1]) / denominator
        
    return corr


def get_dataset_distance_matrix(clustering_dict):
    num_datasets = len(list(clustering_dict.values())[0])
    num_residues = len(clustering_dict)
    
    dataset_connectivity_matrix = np.zeros((num_datasets, num_datasets))
    
    for residue_id, residue_clustering in clustering_dict.items():
        
        for x, cluster_index_x in enumerate(residue_clustering.values()):
            for y, cluster_index_y in enumerate(residue_clustering.values()):
                if cluster_index_x == cluster_index_y:
                    dataset_connectivity_matrix[x,y] = dataset_connectivity_matrix[x,y] + 1
    
    return dataset_connectivity_matrix / num_residues
    
    
def save_num_clusters_bar_plot(clustering_dict, plot_file):
    
    dtag_list = list(list(clustering_dict.values())[0].keys())

    fig, ax = plt.subplots(figsize=(20, 0.1*len(clustering_dict)))
    
    x = np.arange(len(clustering_dict))
    y = [np.unique([cluster_id for cluster_id in cluster_dict.values()]).size for cluster_dict in clustering_dict.values()]
    labels = [dtag.dtag for dtag in dtag_list]
    plt.bar(x, y)
    plt.xticks(x, labels)
    fig.savefig(str(plot_file))
    fig.clear()
    plt.close(fig)

    
def main():
        
    ###################################################################
    # # Configuration
    ###################################################################
    print("Getting config")
    # args = Args.from_cmd(arg_list)
    args = Args(
        Path("/dls/science/groups/i04-1/conor_dev/baz2b_test/data"),
            #    Path("/dls/science/groups/i04-1/conor_dev/experiments/test_local_cluster"),
                Path("/dls/science/groups/i04-1/conor_dev/experiments/test_local_cluster_mask_baz2ba"),
               "*.dimple.pdb" ,
               "*.dimple.mtz",
               )

    # args = Args(Path("/dls/labxchem/data/2017/lb18145-17/processing/analysis/initial_model/"),
    #         # Path("/dls/science/groups/i04-1/conor_dev/experiments/test_local_cluster"),
    #         # Path("/dls/science/groups/i04-1/conor_dev/experiments/test_local_cluster_gauss_diag"),
    #         # Path("/dls/science/groups/i04-1/conor_dev/experiments/test_local_cluster_gauss_diag_no_multi"),
    #         # Path("/dls/science/groups/i04-1/conor_dev/experiments/test_local_cluster_gauss_diag_fast"),
    #         "dimple.pdb" ,
    #         "dimple.mtz",
    #         )

    print("Getting multiprocessor")
    mapper = JoblibMapper.initialise()


    ###################################################################
    # # Pre-pandda
    ###################################################################

    # Get datasets
    print("Loading datasets")
    datasets_initial: Datasets = Datasets.from_data_dirs(args.data_dirs,
                                                        args.pdb_regex,
                                                        args.mtz_regex,
                                                        )
    

    # Initial filters
    print("Filtering invalid datasaets")
    datasets_invalid: Datasets = datasets_initial.remove_invalid_structure_factor_datasets(args.structure_factors)

    datasets_low_res: Datasets = datasets_invalid.remove_low_resolution_datasets(
        args.low_resolution_completeness)

    datasets_rfree: Datasets = datasets_low_res.remove_bad_rfree(args.max_rfree)

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

    min_res = datasets[sorted_dtags[-1]].reflections.resolution()


    print("Truncating datasets")
    truncated_datasets: Datasets = datasets.truncate(resolution=min_res,
                                                                structure_factors=args.structure_factors,
                                                                )



    clustering_dict = {}
    # Iterate over residues
    for residue_id in reference.dataset.structure.protein_residue_ids():
        
        print((
            f"Working on residue: {residue_id}"
        ))

        samples = {dtag:
        sample_residue(
                    truncated_datasets[dtag],
                    grid,
                residue_id,
                    alignments[dtag],
                    args.structure_factors, 
                    args.sample_rate,     
        )
        for dtag in datasets
            }

        correlation_matrix = np.zeros((len(samples), len(samples)))
        
        for x, reference_sample in enumerate(samples.values()):

            reference_sample_copy_unscaled = reference_sample

            reference_sample_copy = reference_sample_copy_unscaled

            reference_sample_mask = np.zeros(reference_sample_copy.shape)
            reference_sample_mask[reference_sample_copy > 1.0] = 1

            for y, sample in enumerate(samples.values()):

                sample_copy = sample.copy()

                sample_mask = np.zeros(sample_copy.shape)
                sample_mask[sample_copy > 1.0] = 1
                
                correlation_matrix[x, y] = get_corr(reference_sample_mask, sample_mask, x==y)

        
        save_correlation_plot(correlation_matrix, 
                              args.out_dir / f"{residue_id}_correlation.png",
                              )
        
        linkage = get_linkage_from_correlation_matrix(correlation_matrix)
        
        cluster_ids = cluster_linkage(linkage, 0.4)
        
        print(cluster_ids)
        
        save_dendrogram_plot(linkage, 
                             labels=[dtag.dtag for dtag in samples.keys()], 
                             dendrogram_plot_file=args.out_dir / f"{residue_id}_dendrogram.png")
        
        clustering_dict[residue_id] = {dtag: cluster_id for dtag, cluster_id in zip(samples.keys(), cluster_ids)}
        
        print(f"Found {np.unique(cluster_ids).size} clusters")


    # #############################
    # Summary
    # #############################
    # # Get dtag list
    dtag_list = list(list(clustering_dict.values())[0].keys())

    # # Connectivity
    dataset_connectivity_matrix = get_dataset_distance_matrix(clustering_dict)
    save_correlation_plot(dataset_connectivity_matrix, 
                            args.out_dir / f"global_connectivity_correlation.png",
                            )
    connectivity_linkaged = get_linkage_from_correlation_matrix(dataset_connectivity_matrix)
    idx = cluster_linkage(connectivity_linkaged, 0.75)
    print(idx)
    print(zip(dtag_list, idx))

    # # Summary plots
    save_dendrogram_plot(connectivity_linkaged, 
                         labels=[dtag.dtag for dtag in dtag_list], 
                         dendrogram_plot_file=args.out_dir / f"global_connectivity_dendrogram.png",
                         )
    
    save_num_clusters_bar_plot(clustering_dict, args.out_dir / f"global_residue_cluster_bar.png")
                
                
if __name__ == "__main__":
    main()