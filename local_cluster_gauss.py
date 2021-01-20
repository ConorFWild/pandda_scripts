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

# from sklearn.decomposition import PCA
# from sklearn.manifold import TSNE
# from hdbscan import HDBSCAN


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

# from scipy import stats
from numpy import random
from sklearn import mixture
from sklearn.neighbors import DistanceMetric



@dataclasses.dataclass()
class Args:
    data_dirs: Path
    out_dir: Path
    pdb_regex: str
    mtz_regex: str
#     structure_factors: StructureFactors = StructureFactors.from_string("2FOFCWT,PH2FOFCWT")
    structure_factors: StructureFactors = StructureFactors.from_string("FWT,PHWT")
    low_resolution_completeness: float = 4.0
    max_rfree: float = 0.4
    max_wilson_plot_z_score: float = 5.0
    max_rmsd_to_reference: float = 1.5
    outer_mask: float = 6.0
    inner_mask_symmetry: float = 3.0
    sample_rate: float = 3.0
    num_components: int = 5
    
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
    
def main():
        
    ###################################################################
    # # Configuration
    ###################################################################
    print("Getting config")
    # args = Args.from_cmd(arg_list)
    # args = Args(Path("/dls/science/groups/i04-1/conor_dev/baz2b_test/data"),
    #            Path("/dls/science/groups/i04-1/conor_dev/experiments/test_local_cluster"),
    #            "*.dimple.pdb" ,
    #            "*.dimple.mtz",
    #            )

    args = Args(Path("/dls/labxchem/data/2017/lb18145-17/processing/analysis/initial_model/"),
            Path("/dls/science/groups/i04-1/conor_dev/experiments/test_local_cluster"),
            "dimple.pdb" ,
            "dimple.mtz",
            )

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

                
    # indexes = np.arange(0, len(datasets_initial), 1)
    # selection = np.random.choice(indexes, size=200, replace=False)
    # keys = list(datasets_initial.datasets.keys())
    # values = list(datasets_initial.datasets.values())
    # new_datasets = {}
    # for index in selection:
    #     key = keys[index]
    #     value = values[index]
    #     new_datasets[key] = value

    # datasets_initial = Datasets(new_datasets)
    

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



    results = {}
    
    # Iterate over residues
    for residue_id in reference.dataset.structure.protein_residue_ids():
        results[residue_id] = {}
        
        print((
            f"Working on residue: {residue_id}"
        ))
        
        selection_list_list = mapper(
            pandda_types.delayed(sample_residue)(
                    truncated_datasets[dtag],
                    grid,
                residue_id,
                    alignments[dtag],
                    args.structure_factors, 
                    args.sample_rate,     
                )
            for dtag
            in truncated_datasets
            )

        dtag_array = np.array([dtag for dtag in truncated_datasets])
            
        sample_by_features = np.vstack(selection_list_list)
        print(sample_by_features.shape)

        # Iterate over different models
        for i in range(args.num_components):
            n_components = i + 1
            results[residue_id][n_components] = {}
            
            # Fit the gmm
            gm = mixture.GaussianMixture(n_components=n_components)
            start_fit_time = time.time()
            gm.fit(sample_by_features)
            finish_fit_time = time.time()
            print(f"Fit gmm in: {finish_fit_time-start_fit_time}")
            
            # Find the aic of the model
            aic = gm.aic(sample_by_features)
            print(f"aic: {aic}")
            
            # Pull out the means and covaraiances
            means = gm.means_
            covs = gm.covariances_
            
            results[residue_id][n_components]["aic"] = aic

            # Get the condifence intervals
            for j, mean, cov in zip(range(len(means)), means, covs):
                results[residue_id][n_components][j] = {}
                
                # Get the distance metric
                dist = DistanceMetric.get_metric('mahalanobis', V=cov)
                
                # Get the distance to known points
                distances = dist.pairwise(sample_by_features, mean.reshape(1,mean.size))
                
                # Get the expected distance to points
                sample = random.multivariate_normal(mean=mean, cov=cov, size = 1000)
                sample_distances = dist.pairwise(sample, mean.reshape(1,mean.size))
                cutoff = np.quantile(sample_distances, 0.9)

                print(f"Component {j}")
                print(f"Max distance: {np.max(distances)}")
                print(f"Cutoff: {cutoff}")
                print(f"Number of distances less than {cutoff}: {distances[distances<cutoff].size}")

                # Create a plot
                fig = ff.create_distplot([distances.flatten()], ["distances"], bin_size=0.5)
                fig.write_image(str(args.out_dir / f"test_{residue_id}_{n_components}_{j}_dist.png"),
                                engine="kaleido", 
                                width=2000,
                                height=1000,
                                scale=1)  
                
                results[residue_id][n_components][j]["num_distances"] = distances[distances<cutoff].size
                results[residue_id][n_components][j]["dtags"] = [dtag.dtag for dtag in dtag_array[(distances<cutoff).flatten()]]
                
if __name__ == "__main__":
    main()