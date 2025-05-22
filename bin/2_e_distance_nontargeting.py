#!/usr/bin/env python
# coding: utf-8

import sys
import os

import pandas as pd
import numpy as np
import scanpy as sc

from sklearn.cluster import KMeans
import scipy.stats
from itertools import combinations
from sklearn.metrics import pairwise_distances

import seaborn as sns
import matplotlib.pyplot as plt
plt.rcParams.update({'axes.labelsize' : 'large',
                     'pdf.fonttype' : 42
                    })
from matplotlib.backends.backend_pdf import PdfPages

from tqdm import tqdm
from tqdm.contrib.concurrent import process_map

import gc
import warnings
import time
import pickle
import json

from sklearn.metrics import pairwise_distances
from multiprocessing import Pool
import torch

from importlib import reload
import util_functions
import energy_distance_calc



json_fp = sys.argv[1]
with open(json_fp, 'r') as fp:
    config = json.load(fp)

device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

(pca_df,gRNA_dict) = util_functions.load_files(config["input_data"]["h5ad_file"],
                                               config["input_data"]["sgRNA_file"],
                                               os.path.join(config["output_file_name_list"]["OUTPUT_FOLDER"],
                                                            config["output_file_name_list"]["pca_table"]),
                                               os.path.join(config["output_file_name_list"]["OUTPUT_FOLDER"],
                                                            config["output_file_name_list"]["gRNA_dict"]),
                                               overwrite=False
                                              )

sgRNA_outlier_df = pd.read_csv(os.path.join(config["output_file_name_list"]["OUTPUT_FOLDER"],
                                                  config["output_file_name_list"]["targeting_outlier_table"]),
                                     index_col=0)

nontargeting_outlier_df = pd.read_csv(os.path.join(config["output_file_name_list"]["OUTPUT_FOLDER"],
                                                         config["output_file_name_list"]["non_targeting_outlier_table"]),
                                            index_col=0)


clear_sgRNA_list = sgRNA_outlier_df[sgRNA_outlier_df["pval_outlier"]>0.05].index.tolist()
clear_nt_sgRNA_list = nontargeting_outlier_df[nontargeting_outlier_df["pval_outlier"]>0.05].index.tolist()


annotation_df = pd.read_csv(os.path.join(config["output_file_name_list"]["OUTPUT_FOLDER"],
                                         config["output_file_name_list"]["annotation_file"]),index_col=0)

gRNA_region_dict = util_functions.get_gRNA_region_dict(annotation_df,gRNA_dict)

gRNA_region_clear_dict = {}

for key in gRNA_region_dict.keys():
    gRNA_list_tmp = [x for x in gRNA_region_dict[key] if x in clear_sgRNA_list]
    if len(gRNA_list_tmp)!=0:
        gRNA_region_clear_dict[key] = [x for x in gRNA_region_dict[key] if x in clear_sgRNA_list]


# In[13]:


cell_per_region_dict = {}
for key in gRNA_region_clear_dict.keys():
    cell_list_tmp = [gRNA_dict[i] for i in gRNA_region_clear_dict[key]]
    cell_list_tmp = np.concatenate(cell_list_tmp)
    cell_per_region_dict[key] = np.unique(cell_list_tmp)



print("Total:",len(gRNA_region_clear_dict.keys()))


permute_per_bg = config["permutation_test"]["permute_per_bg"]
num_of_bg = config["permutation_test"]["num_of_bg"]
non_target_pick = config["permutation_test"]["non_target_pick"]
batch_num_basic = config["permutation_test"]["batch_num_basic"]

d_tmp = np.log10(0.1/(permute_per_bg*num_of_bg))
d_tmp = np.round(d_tmp)
d_tmp = np.power(10,d_tmp)


non_target_cell_name = [gRNA_dict[i] for i in clear_nt_sgRNA_list]

#concatenate all list together
non_target_cell_name = sum(non_target_cell_name,[])

#unique cell id list
non_target_cell_name = sorted(list(set(non_target_cell_name)))


#fix random cell with seed
np.random.seed(0)

non_target_cell_name_reduced = []

for i in range(num_of_bg):
    non_target_cell_name_reduced.append(np.random.choice(non_target_cell_name,non_target_pick))


output_file = os.path.join(config["output_file_name_list"]["OUTPUT_FOLDER"],
                           config["output_file_name_list"]["edist_pvalue_table"])


## Energy distance calculation

### Helper functions
def get_type_target(target_name,annotation_df,annotation_col="intended_target_name"):
    return annotation_df.loc[(annotation_df[annotation_col]==target_name),"type"].values[0]

def adjust_processing_parameters(target_cell_names, default_batch_num):
    """
    Adjusts batch number and performs downsampling based on the number of target cells.

    Args:
        target_cell_names (list or np.ndarray): List of target cell names.
        default_batch_num (int): The default batch number to use.

    Returns:
        tuple: A tuple containing:
            - adjusted_cell_names (list or np.ndarray): Original or downsampled cell names.
            - batch_num (int): The adjusted batch number.
    """
    cell_count = len(target_cell_names)
    adjusted_cell_names = target_cell_names
    processing_flag = 0 # Default flag

    if cell_count > 5000:
        batch_num = 40
        # Downsample to 5000 cells without replacement
        adjusted_cell_names = np.random.choice(target_cell_names, 5000, replace=False)
        print(f"Target region with initially {cell_count} cells is downsampled to 5000.")
    elif cell_count > 2000:
        batch_num = 50
    elif cell_count > 1000:
        batch_num = 100
    else:
        batch_num = default_batch_num # Use default if below thresholds

    # Note: The 'processing_flag' seems unused later in the provided snippet,
    # but is retained to match the original logic's variable assignments.
    return adjusted_cell_names, batch_num

def remove_overlap_cells(target_cells, non_target_cells):
    """
    Removes overlapping cell names between target and non-target sets.

    Args:
        target_cells (list or set): Collection of target cell names.
        non_target_cells (list or set): Collection of non-target cell names.

    Returns:
        tuple: A tuple containing:
            - target_cells_clear (list): Target cells with overlap removed.
            - non_target_cells_clear (list): Non-target cells with overlap removed.
    """
    target_set = set(target_cells)
    non_target_set = set(non_target_cells)

    overlap = target_set & non_target_set
    target_cells_clear = list(target_set - overlap)
    non_target_cells_clear = list(non_target_set - overlap)

    return target_cells_clear, non_target_cells_clear

def run_permutation_test_with_fallback(pbar, pca_df, target_cells, non_target_cells,
                                     gpu_device, batch_num, permutations,
                                     current_iter_info, total_cell_info):
    """
    Runs the energy distance permutation test, attempting GPU first and falling back to CPU.

    Args:
        pbar (tqdm): The progress bar instance to update postfix.
        pca_df (pd.DataFrame): DataFrame containing PCA results.
        target_cells (list): List of target cell names (cleared of overlap).
        non_target_cells (list): List of non-target cell names (cleared of overlap).
        gpu_device (torch.device or str): The primary device (GPU) to try.
        batch_num (int): Batch number for calculation.
        permutations (int): Number of permutations for the test.
        current_iter_info (any): Information about the current iteration for postfix.
        total_cell_info (any): Information about the total cells for postfix.


    Returns:
        tuple: A tuple containing:
            - obs_edist (float or None): Observed energy distance, or None if calculation failed.
            - e_dist_list (list or None): List of energy distances from permutations, or None.
    """
    obs_edist, e_dist_list = None, None
    mode = gpu_device # Assume GPU initially

    try:
        pbar.set_postfix({
            "total cell num": total_cell_info,
            "current iter": current_iter_info,
            "mode": mode
        })
        # Attempt GPU calculation
        obs_edist, e_dist_list = energy_distance_calc.permutation_test(
            pca_df, target_cells, non_target_cells,
            gpu_device, batch_num, permutations
        )
        # print(f"Calculation successful on {mode}") # Optional: uncomment for verbose output
        return obs_edist, e_dist_list

    except Exception as e_gpu:
        print(f"GPU calculation failed: {e_gpu}. Attempting CPU fallback...")
        # Clean up GPU memory before CPU attempt
        gc.collect()
        if torch.cuda.is_available():
             torch.cuda.empty_cache()

        mode = "CPU" # Switch mode for postfix and calculation
        pbar.set_postfix({
            "total cell num": total_cell_info,
            "current iter": current_iter_info,
            "mode": mode
        })

        try:
            # Attempt CPU calculation
            obs_edist, e_dist_list = energy_distance_calc.permutation_test(
                pca_df, target_cells, non_target_cells,
                "cpu", batch_num, permutations # Explicitly use "cpu"
            )
            # print(f"Calculation successful on {mode}") # Optional: uncomment for verbose output
            return obs_edist, e_dist_list

        except Exception as e_cpu:
            # Both GPU and CPU attempts failed
            print(f"CPU calculation also failed: {e_cpu}")
            print("Skipping energy distance calculation for this iteration.")
            # Return None to indicate failure
            return None, None


def save_results(data_dict, output_filename,pval_d=0.00001):
    """
    Saves the results dictionary to a CSV file.

    Args:
        data_dict (dict): The dictionary containing results.
        output_filename (str): Path to the output CSV file.
    """

    for key in data_dict.keys():
        distance_sum = 0
        pval_sum = 0

        for bg_index in range(num_of_bg):
            distance_sum += data_dict[key]["distance_"+str(bg_index)]
            pval_sum += data_dict[key]["pval_"+str(bg_index)]

        data_dict[key]["distance_mean"] = distance_sum/num_of_bg
        data_dict[key]["pval_mean"] = pval_sum/num_of_bg

        data_dict[key]["pval_mean_log"] = -np.log10(data_dict[key]["pval_mean"]+d_tmp)
        data_dict[key]["distance_mean_log"] = np.log10(data_dict[key]["distance_mean"])

    try:
        pd.DataFrame(data_dict).T.to_csv(output_filename)
        # print(f"Results saved to {output_filename}") # Optional: uncomment for confirmation
    except Exception as e:
        print(f"Error saving results to {output_filename}: {e}")

# Initialize results dictionary
pval_dict = {}

# Initialize progress bar
target_regions = list(gRNA_region_dict.keys())
pbar = tqdm(enumerate(target_regions), total=len(target_regions), desc="Processing Regions")

### Main loop iterating through target regions
for target_index, target_name in pbar:
    pbar.set_description(f"Processing {target_name}")
    original_target_cell_names = cell_per_region_dict.get(target_name, []) # Get cells for the region

    if len(original_target_cell_names)==0:
        print(f"Warning: No cells found for target region {target_name}. Skipping.")
        continue

    # Adjust batch size and downsample if necessary
    target_cell_names, batch_num = adjust_processing_parameters(
        original_target_cell_names, batch_num_basic
    )

    # Inner loop for background comparisons
    for bg_index in range(num_of_bg):
        current_bg_non_target_cells = non_target_cell_name_reduced[bg_index]

        # Remove overlapping cells between target and current background
        target_cell_name_clear, non_target_cell_name_clear = remove_overlap_cells(
            target_cell_names, current_bg_non_target_cells
        )

        # Check if lists are empty after removing overlap (can happen if sets are identical)
        if not target_cell_name_clear or not non_target_cell_name_clear:
            print(f"Warning: No unique cells left for comparison between {target_name} and background {bg_index+1}. Skipping.")
            continue

        # Run permutation test with GPU preference and CPU fallback
        obs_edist, e_dist_list = run_permutation_test_with_fallback(
            pbar=pbar,
            pca_df=pca_df,
            target_cells=target_cell_name_clear,
            non_target_cells=non_target_cell_name_clear,
            gpu_device=device,
            batch_num=batch_num,
            permutations=permute_per_bg,
            current_iter_info=f"{bg_index+1}/{num_of_bg}",
            total_cell_info=len(target_cell_names) # Use count after potential downsampling
        )

        # Process results only if calculation was successful
        if obs_edist is not None and e_dist_list is not None:
            # Initialize dict for target_name if first time
            if target_name not in pval_dict:
                pval_dict[target_name] = {}
                pval_dict[target_name]["cell_count"] = len(target_cell_names)
                pval_dict[target_name]["type"] = get_type_target(target_name,annotation_df)

            # Store observed distance and calculate p-value
            pval_dict[target_name][f"distance_{bg_index}"] = obs_edist.item() # Use .item() for scalar tensor
            # Calculate p-value: proportion of permutation distances >= observed distance
            p_value = (e_dist_list >= obs_edist).sum().item() / permute_per_bg
            pval_dict[target_name][f"pval_{bg_index}"] = p_value
        else:
            # Handle failure case (optional: store NaN or skip)
            if target_name not in pval_dict:
                pval_dict[target_name] = {}
            pval_dict[target_name][f"distance_{bg_index}"] = np.nan # Store NaN if failed
            pval_dict[target_name][f"pval_{bg_index}"] = np.nan      # Store NaN if failed
            print(f"Stored NaN for {target_name}, background {bg_index} due to calculation failure.")

    # Save results periodically (e.g., every 10 regions)
    if (target_index + 1) % 10 == 0: # Check using target_index + 1 if saving after processing 10, 20,...
        print(f"\nSaving intermediate results after processing {target_index + 1} regions...")
        save_results(pval_dict, output_file, d_tmp)
        print("Intermediate results saved.\n")


# Final save after processing all regions
print("\nProcessing complete. Saving final results...")
save_results(pval_dict, output_file, d_tmp)
print(f"Final results saved to {output_file}")
