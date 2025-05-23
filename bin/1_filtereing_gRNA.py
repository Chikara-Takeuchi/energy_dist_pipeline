#!/usr/bin/env python
# coding: utf-8



import sys
import os

import pandas as pd
import numpy as np
import scanpy as sc

from sklearn.cluster import KMeans
from scipy.stats import hypergeom
from sklearn.metrics import pairwise_distances
from itertools import combinations


import seaborn as sns
import matplotlib.pyplot as plt
plt.rcParams.update({'axes.labelsize' : 'large',
                     'pdf.fonttype':42
                    })
from matplotlib.backends.backend_pdf import PdfPages

from tqdm import tqdm
from tqdm.contrib.concurrent import process_map

import multiprocessing as mp

import gc
import warnings
import time
import pickle
import json
import math

import torch

from importlib import reload
import util_functions
import energy_distance_calc


### Load PCA matric and gRNA matrix, and related information
json_fp = sys.argv[1]
with open(json_fp, 'r') as fp:
    config = json.load(fp)

device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

annotation_df = pd.read_csv(config["input_data"]["annotation_file"],index_col=0)


(pca_df,gRNA_dict) = util_functions.load_files(config["input_data"]["h5ad_file"],
                                               config["input_data"]["sgRNA_file"],
                                               os.path.join(config["output_file_name_list"]["OUTPUT_FOLDER"],
                                                            config["output_file_name_list"]["pca_table"]),
                                               os.path.join(config["output_file_name_list"]["OUTPUT_FOLDER"],
                                                            config["output_file_name_list"]["gRNA_dict"]),
                                               overwrite=config["output_file_name_list"]["OVERWRITE_PCA_DICT"]
                                              )


gRNA_region_dict = {}
count_region_dict = {}

for index,row in annotation_df.iterrows():
    if row.guide_id in gRNA_dict.keys():
        if row.intended_target_name in gRNA_region_dict.keys():
            gRNA_region_dict[row.intended_target_name] += [row.guide_id]
        else:
            gRNA_region_dict[row.intended_target_name] = [row.guide_id]

for key in gRNA_region_dict.keys():
    gRNA_region_dict[key] = np.unique(gRNA_region_dict[key])
    count_region_dict[key] = len(gRNA_region_dict[key])


### Calculate e-distance between gRNA per target and non-targeting

# Helper functions
def generate_sgrna_group_combinations(gRNA_list_target, combi_count):
    """
    Generates unique pairs of sgRNA groups based on the input list.

    Handles two cases based on the total number of gRNAs:
    - If > 6: Uses a fixed combination size ('combi_count').
    - If <= 6: Generates combinations of varying sizes.

    Args:
        gRNA_list_target (np.array): Array of sgRNA names for the current target region.
        combi_count (int): The number of items to choose when total gRNAs > 6.

    Returns:
        list: A list of unique tuples, where each tuple contains two lists
              representing a pair of sgRNA groups. e.g., [(['gA1'], ['gA2', 'gA3']), ...]
    """
    total_combis_tmp = []
    total_gRNA_num = len(gRNA_list_target)
    gRNA_list_target_tuple = tuple(gRNA_list_target) # Use tuple for combinations

    if total_gRNA_num > 6:
        # Case 1: More than 6 gRNAs, use fixed combination size 'combi_count'
        if total_gRNA_num < combi_count:
            print(f"Warning: total_gRNA_num ({total_gRNA_num}) < combi_count ({combi_count}). Skipping combination generation for this target.")
            return []
        # Generate combinations of size 'combi_count'
        for combis_tmp_tuple in combinations(gRNA_list_target_tuple, combi_count):
            combis_tmp_set = set(combis_tmp_tuple)
            # Split these 'combi_count' gRNAs into two non-empty groups
            for first_group_num in range(1, combi_count):
                for first_group_tuple in combinations(combis_tmp_tuple, first_group_num):
                    # Efficiently find the second group using set difference
                    second_group_set = combis_tmp_set - set(first_group_tuple)
                    if second_group_set: # Ensure the second group is not empty
                        total_combis_tmp.append((list(first_group_tuple), list(second_group_set)))

    else:
        # Case 2: 6 or fewer gRNAs, generate all pairs of non-empty, disjoint subsets
        for num_count_1 in range(1, total_gRNA_num): # Size of the first group
            for group1_tuple in combinations(gRNA_list_target_tuple, num_count_1):
                # Find remaining gRNAs efficiently
                remaining_gRNAs = tuple(g for g in gRNA_list_target_tuple if g not in group1_tuple)
                if not remaining_gRNAs: # Need remaining gRNAs to form the second group
                    continue
                # Generate the second group from the remaining gRNAs
                for num_count_2 in range(1, len(remaining_gRNAs) + 1): # Size of the second group
                     for group2_tuple in combinations(remaining_gRNAs, num_count_2):
                            total_combis_tmp.append((list(group1_tuple), list(group2_tuple)))

    # Ensure uniqueness of pairs (consider order within pair irrelevant for uniqueness check initially)
    # Convert lists to tuples for hashing, sort internal lists first, then the pair
    # E.g. (['B', 'A'], ['C']) becomes (('A', 'B'), ('C',))
    # E.g. (['C'], ['A', 'B']) becomes (('A', 'B'), ('C',))
    seen_combinations = set()
    unique_combis = []
    for group1, group2 in total_combis_tmp:
        tuple1 = tuple(sorted(group1))
        tuple2 = tuple(sorted(group2))
        # Sort the pair of tuples to handle cases like (g1, g2) vs (g2, g1) being the same split
        frozen_pair = tuple(sorted((tuple1, tuple2)))
        if frozen_pair not in seen_combinations:
            seen_combinations.add(frozen_pair)
            # Add back in the original list format
            unique_combis.append((group1, group2))

    # Note: The original code called util_functions.get_unique_list AFTER this.
    # If get_unique_list does something more complex than the uniqueness logic here,
    # you might need to adjust or call it after this function.
    # Assuming the logic above correctly finds unique pairs.
    return unique_combis


def get_cells_for_sgrna_groups(group1, group2, gRNA_dict):
    """
    Retrieves and combines unique cell identifiers for two groups of sgRNAs.

    Args:
        group1 (list): List of sgRNA names in the first group.
        group2 (list): List of sgRNA names in the second group.
        gRNA_dict (dict): Dictionary mapping sgRNA names to lists/arrays of cell identifiers.

    Returns:
        tuple: (np.array, np.array) containing unique cell identifiers for group1 and group2.
               Returns (None, None) if a gRNA name is not found in gRNA_dict.
    """
    try:
        cells1 = np.unique(np.concatenate([gRNA_dict[name] for name in group1])) if group1 else np.array([])
        cells2 = np.unique(np.concatenate([gRNA_dict[name] for name in group2])) if group2 else np.array([])
        return cells1, cells2
    except KeyError as e:
        print(f"Error: sgRNA name {e} not found in gRNA_dict.")
        return None, None
    except Exception as e:
        print(f"An error occurred retrieving cells: {e}")
        return None, None

def calculate_distance(X, cell_test1, cell_test2, device):
    """
    Calculates the distance between two cell groups using permutation_test.
    Attempts GPU calculation first, falls back to CPU if memory error occurs.

    Args:
        X (np.array or similar): The data matrix (e.g., PCA coordinates).
        cell_test1 (np.array): Array of cell identifiers for the first group.
        cell_test2 (np.array): Array of cell identifiers for the second group.
        device (str): The primary device to use ('cuda', 'cpu', etc.).

    Returns:
        float: The calculated distance, or -1.0 if calculation fails.
    """
    if cell_test1 is None or cell_test2 is None: # Check if cell retrieval failed
        return -1.0
    if len(cell_test1) == 0 or len(cell_test2) == 0:
        # Handle cases where one group has no cells (might happen if gRNA mapping is empty)
        # Distance is ill-defined or could be considered infinite/maximal. Return -1 as indicator.
        # print(f"Warning: Empty cell group for combination {index_combi}. Skipping distance calculation.")
        return -1.0

    mode = "GPU" if "cuda" in str(device) else "CPU" # Initial mode assumption
    try:
        # Assuming permutation_test takes numpy arrays for cell indices/names mapped to X
        # Ensure X is indexed correctly based on cell_test1/cell_test2 content
        # The original call used cell names/indices directly, assuming X is indexed accordingly.
        obs_edist = energy_distance_calc.permutation_test(X, cell_test1, cell_test2, device,
                                                     1, 1, return_permute=False).cpu()
        return obs_edist.item()
    except Exception as e_gpu: # Catch specific OOM error if possible, otherwise generic Exception
        # Check if it looks like an Out Of Memory error (heuristic)
        if 'memory' in str(e_gpu).lower():
            mode = "OOM -> CPU"
            try:
                # Fallback to CPU
                obs_edist = energy_distance_calc.permutation_test(X, cell_test1, cell_test2, "cpu",
                                                             1, 1, return_permute=False).cpu()
                return obs_edist.item()
            except Exception as e_cpu:
                # Handle CPU failure (e.g., data still too large even for CPU RAM)
                print(f"\nError: CPU calculation failed for combi {index_combi+1} after GPU OOM. Data too large? Error: {e_cpu}")
                return -1.0
        else:
            # Handle other non-OOM GPU errors
            print(f"\nError: GPU calculation failed for combi {index_combi+1} (not OOM). Error: {e_gpu}")
            return -1.0


def format_results_dataframe(res, total_combis, gRNA_list_target):
    """
    Creates a DataFrame from the distance results, sorts it, adds ranks,
    and includes boolean flags for sgRNA membership in each combination.

    Args:
        res (list): List of calculated distances (-1.0 indicates failure).
        total_combis (list): List of tuples, each containing two lists (sgRNA groups).
        gRNA_list_target (np.array): Array of all sgRNA names for the current target.

    Returns:
        pd.DataFrame: A formatted DataFrame with results.
    """
    if len(res) != len(total_combis):
        print(f"Warning: Mismatch between results count ({len(res)}) and combinations count ({len(total_combis)}). DataFrame might be incorrect.")
        # Attempt to proceed with the shorter length
        min_len = min(len(res), len(total_combis))
        res = res[:min_len]
        total_combis = total_combis[:min_len]

    result_df = pd.DataFrame({
        "e_dist": res,
        "combis": total_combis
    })

    # Filter out failed calculations if needed (optional, depends on desired output)
    # result_df = result_df[result_df["e_dist"] != -1.0]

    result_df = result_df.sort_values(by="e_dist").reset_index(drop=True)
    result_df["rank"] = result_df.index

    # Add boolean columns for each sgRNA
    for gRNA_name_tmp in gRNA_list_target:
        result_df[gRNA_name_tmp] = result_df["combis"].apply(
            lambda x: (gRNA_name_tmp in x[0]) or (gRNA_name_tmp in x[1]) if isinstance(x, (list, tuple)) and len(x) == 2 else False
        )

    return result_df


# In[12]:


# --- Main Processing Loop ---

# gRNA_region_dict: Dict mapping region name to list/array of sgRNAs
# gRNA_dict: Dict mapping sgRNA name to list/array of cells
# X: Data matrix (e.g., PCA) indexed appropriately
# device: Computation device ('cuda', 'cpu')
# combi_count: Parameter for combination generation (e.g., 6)
# result_df_dict: Dictionary to store results per target
combi_count = 6

result_df_dict = {} # Initialize result dictionary

test_region_np = np.array(list(gRNA_region_dict.keys()))
test_region_np = test_region_np[~(test_region_np=="non-targeting")]

print(f"Processing {len(test_region_np)} target regions...")

# Outer loop iterates through target regions
pbar_targets = tqdm(test_region_np, desc="Overall Target Regions")
for target in pbar_targets:
    pbar_targets.set_postfix({
            "Current Target": target
        })

    gRNA_list_target = np.array(gRNA_region_dict[target])
    if len(gRNA_list_target) < 2:
        print(f"Skipping target '{target}': Needs at least 2 gRNAs, found {len(gRNA_list_target)}.")
        continue

    # 1. Generate combinations of sgRNA groups
    # print(f"\nGenerating combinations for target: {target}") # Optional verbose print
    total_combis = generate_sgrna_group_combinations(gRNA_list_target, combi_count)

    if not total_combis:
        print(f"Skipping target '{target}': No valid combinations generated.")
        continue

    # print(f"Generated {len(total_combis)} unique combinations for {target}.") # Optional verbose print

    # 2. Calculate distance for each combination
    res = []
    # Use a nested tqdm for combinations within a target, leaving the outer one clean
    #pbar_combis = tqdm(enumerate(total_combis), total=len(total_combis), desc=f"Target '{target}' Combinations", leave=False)
    for index_combi, (combi_test1, combi_test2) in enumerate(total_combis):
        # Get cells for the current combination
        cell_test1, cell_test2 = get_cells_for_sgrna_groups(combi_test1, combi_test2, gRNA_dict)

        # Try calculating on the primary device (potentially GPU)
        # Calculate distance with GPU/CPU fallback
        distance = calculate_distance(pca_df, cell_test1, cell_test2, device)
        res.append(distance)

    # 3. Format results into a DataFrame
    # print(f"\nFormatting results for target: {target}") # Optional verbose print
    result_df = format_results_dataframe(res, total_combis, gRNA_list_target)

    # Store the result DataFrame
    result_df_dict[target] = result_df
    # print(f"Finished processing target: {target}") # Optional verbose print

print("\nAll target regions processed.")
# Now result_df_dict contains the DataFrames for each target region.


target_arr = np.array(list(result_df_dict.keys()))


def determine_batch_size(total_cell_num, batch_num_basic):
    """
    Determines the appropriate batch size based on the total number of cells.

    Args:
        total_cell_num (int): The total number of cells for the current target.
        batch_num_basic (int): The default base batch number.

    Returns:
        int: The calculated batch size (at least 1).
    """
    if total_cell_num > 20000:
        batch_num = 1 # Smallest batch size for very large datasets
    elif total_cell_num > 5000:
        batch_num = batch_num_basic // 4 # Reduced batch size for large datasets
    elif total_cell_num > 300:
        batch_num = batch_num_basic // 2 # Reduced batch size for medium datasets
    else:
        batch_num = batch_num_basic # Basic batch size for smaller datasets
    # Ensure batch size is at least 1
    return max(1, batch_num)

def run_disco_test(X, total_cell_list, device, batch_num, total_permute_disco, target_name):
    """
    Runs the disco test, attempting GPU first and falling back to CPU.

    Args:
        X (array-like): The data matrix (e.g., PCA coordinates).
        total_cell_list (list): A list where each element is a list/array of cell IDs
                                corresponding to an sgRNA group.
        device (str): The primary computation device ('cuda', 'cpu', etc.).
        batch_num (int): The batch number for the test.
        total_permute_disco (int): The total number of permutations used in the disco test.
        target_name (str): The name of the target region (for logging).

    Returns:
        float: The calculated disco test p-value, or np.nan if the test fails completely.
    """
    # Removed print statement for batch size here, logged before calling
    try:
        # Attempt calculation on the primary device (potentially GPU)
        # print(f"    Attempting disco test on device: {device} for {target_name}", flush=True) # Optional verbose
        obs_fvalue, fvalue_list = energy_distance_calc.disco_test(X, total_cell_list, device, batch_num=batch_num)
        obs_fvalue = obs_fvalue.numpy(),
        fvalue_list = fvalue_list.numpy()

        disco_pvalue = np.sum(fvalue_list > obs_fvalue) / total_permute_disco
        #print(f"    Disco test completed on {device} for {target_name}. p-value: {disco_pvalue:.4f}", flush=True)
        return disco_pvalue
    except Exception as e_gpu:
        error_str = str(e_gpu).lower()
        if 'memory' in error_str or 'cuda' in error_str:
            print(f"    GPU execution failed for {target_name} (likely OOM). Falling back to CPU. Error: {e_gpu}", flush=True)
            try:
                # print(f"    Attempting disco test on CPU for {target_name}", flush=True) # Optional verbose
                obs_fvalue, fvalue_list = util_functions.disco_test(X, total_cell_list, "cpu", batch_num=batch_num)
                disco_pvalue = np.sum(np.asarray(fvalue_list) > obs_fvalue) / total_permute_disco
                #print(f"    Disco test completed on CPU for {target_name}. p-value: {disco_pvalue:.4f}", flush=True)
                return disco_pvalue
            except Exception as e_cpu:
                print(f"    CPU execution also failed for {target_name}. Skipping disco test. Error: {e_cpu}", flush=True)
                return np.nan
        else:
            #print(f"    GPU execution failed for {target_name} (non-memory error). Skipping disco test. Error: {e_gpu}", flush=True)
            return np.nan


def calculate_hypergeometric_pvalue(result_df, gRNA_name, significance_fraction):
    """
    Calculates the p-value for a single sgRNA using the hypergeometric test.
    Tests if the sgRNA is overrepresented in the highest e-dist ranks.

    Args:
        result_df (pd.DataFrame): DataFrame containing sorted e-distance results
                                  and boolean flags for sgRNA membership per combination.
        gRNA_name (str): The name of the sgRNA to calculate the p-value for.
        significance_fraction (float): The fraction of top ranks to consider significant.

    Returns:
        float: The calculated hypergeometric p-value (right-tailed). Returns 1.0
               if input parameters are invalid for the test.
    """
    M = result_df.shape[0] # Population size (total combinations)
    if M == 0:
        # This case should ideally be caught before calling this function
        return 1.0

    n = result_df[gRNA_name].sum() # Successes in population (combinations with the gRNA)
    # Ensure num_sig_diff calculation avoids zero, minimum 1 draw needed for test
    num_sig_diff = max(1, int(M * significance_fraction)) # Number of draws (top ranked combinations)

    # Get top ranks (highest e_dist -> bottom rows if sorted ascending)
    result_df_right = result_df.iloc[-num_sig_diff:]
    x = result_df_right[gRNA_name].sum() # Successes in draw (gRNA presence in top ranks)
    N = num_sig_diff # Number of draws

    # Validate parameters before passing to hypergeom.sf
    if not (0 <= x <= min(n, N) and 0 < N <= M and 0 <= n <= M):
        print(f"    Warning: Invalid parameters for hypergeometric test for {gRNA_name}. "
              f"[M={M}, n={n}, N={N}, x={x}]. Returning p=1.0")
        return 1.0

    # P(X >= x) = 1 - P(X <= x-1) = sf(x-1)
    p_val = hypergeom.sf(x - 1, M, n, N)
    return round(p_val, 4)


# --- Constants ---
# Define constants near the top for clarity and easy modification
DISCO_P_VALUE_THRESHOLD = 0.05
HYPERGEOM_SIGNIFICANCE_FRACTION = 0.5 # Top 10% for hypergeometric test

# --- Main Processing Logic ---

# Assume these variables are defined from previous steps:
# result_df_dict: Dict mapping target name to results DataFrame
# gRNA_region_dict: Dict mapping target region name to list/array of sgRNAs
# gRNA_dict: Dict mapping sgRNA name to list/array of cells
# X: Data matrix (e.g., PCA)
# device: Computation device ('cuda', 'cpu')
# batch_num_basic: Default batch number (integer)
# total_permute_disco: Total permutations for disco test (integer)


# Initialize dictionaries
p_val_dict = {}
disco_val_dict = {}

print("Starting Step 2: Disco Tests and Hypergeometric Outlier Calculation...")

target_keys = list(result_df_dict.keys())
num_targets = len(target_keys)

for i, target in tqdm(enumerate(target_keys),total=len(target_keys)):

    # --- 1. Data Preparation ---
    if target not in result_df_dict or target not in gRNA_region_dict:
        print(f"  Warning: Missing data for target '{target}'. Skipping.")
        continue

    result_df = result_df_dict[target]
    gRNA_list_target = gRNA_region_dict[target] # Renamed for clarity

    if not isinstance(gRNA_list_target, (list, np.ndarray)) or len(gRNA_list_target) == 0:
        print(f"  Warning: Invalid or empty gRNA list for target '{target}'. Skipping.")
        continue

    # Prepare cell list and count total cells, handle errors
    total_cell_list = []
    total_cell_num = 0
    valid_groups = 0 # Count groups that actually have cells
    try:
        cell_lists_temp = []
        for gRNA_name in gRNA_list_target:
            cells = gRNA_dict[gRNA_name] # Raises KeyError if missing
            cell_lists_temp.append(cells)
            if len(cells) > 0:
                 valid_groups += 1
        # Concatenate only non-empty lists if needed, but disco_test might handle list of lists directly
        total_cell_list = cell_lists_temp # Pass the full list structure
        # Calculate total cells more carefully if concatenation is memory intensive
        total_cell_num = sum(len(cells) for cells in total_cell_list)

    except KeyError as e:
        print(f"  Error: sgRNA name {e} not found in gRNA_dict for target '{target}'. Skipping.")
        # Clean up potentially large objects for this iteration before continuing
        del result_df, gRNA_list_target
        gc.collect()
        continue

    #print(f"  Target: {target}, Num gRNAs: {len(gRNA_list_target)}, Total Cells: {total_cell_num}, Groups with Cells: {valid_groups}")

    # --- 2. Determine Batch Size ---
    batch_num = determine_batch_size(total_cell_num, config["gRNA_filtering"]["batch_num_basic"])
    #print(f"  Determined batch size: {batch_num}")

    # --- 3. Run Disco Test ---
    if valid_groups < 2:
        print(f"  Skipping disco test for '{target}': Fewer than 2 sgRNAs have associated cells.")
        disco_pvalue = np.nan
    else:
        # Assuming disco_test can handle the list of cell lists directly
        disco_pvalue = run_disco_test(pca_df, total_cell_list, device, batch_num, config["gRNA_filtering"]["total_permute_disco"], target)

    disco_val_dict[target] = disco_pvalue

    # --- 4. Calculate Individual sgRNA p-values (Conditional) ---
    # Check disco test outcome
    if math.isnan(disco_pvalue) or disco_pvalue > DISCO_P_VALUE_THRESHOLD or len(gRNA_list_target) == 2:
        if math.isnan(disco_pvalue):
             log_msg = f"Disco test failed or skipped for '{target}'."
        elif disco_pvalue > DISCO_P_VALUE_THRESHOLD:
             log_msg = f"Disco test p-value ({disco_pvalue:.4f}) > {DISCO_P_VALUE_THRESHOLD} for '{target}'."
        else: # len == 2 case
             log_msg = f"Target '{target}' has only 2 sgRNAs."
        #print(f"  {log_msg} Assigning default p-value (1.0) to its sgRNAs.")

        for gRNA_name_tmp in gRNA_list_target:
            p_val_dict[gRNA_name_tmp] = 1.0

    else:
        # Disco test significant, proceed with hypergeometric tests
        #print(f"  Disco test significant (p={disco_pvalue:.4f}). Calculating hypergeometric p-values for sgRNAs in '{target}'.")
        if result_df.empty:
            print(f"  Warning: result_df for target '{target}' is empty. Cannot calculate hypergeometric p-values. Assigning 1.0.")
            for gRNA_name_tmp in gRNA_list_target:
                p_val_dict[gRNA_name_tmp] = 1.0
        else:
            # Check if result_df contains expected columns
            missing_cols = [g for g in gRNA_list_target if g not in result_df.columns]
            if missing_cols:
                print(f"  Error: Missing gRNA columns in result_df for target '{target}': {missing_cols}. Assigning 1.0 p-value.")
                for gRNA_name_tmp in gRNA_list_target:
                    p_val_dict[gRNA_name_tmp] = 1.0
            else:
                # Calculate p-value for each sgRNA
                for gRNA_name_tmp in gRNA_list_target:
                    p_value = calculate_hypergeometric_pvalue(result_df, gRNA_name_tmp, HYPERGEOM_SIGNIFICANCE_FRACTION)
                    p_val_dict[gRNA_name_tmp] = p_value
                    # print(f"    p-value for {gRNA_name_tmp}: {p_value}") # Optional verbose print
    # gc.collect() # Uncomment if memory issues are severe, but adds overhead

# --- 5. Final Output ---
print("\nFinished processing all targets.")
print("Aggregating results and saving to CSV...")

outlier_df = pd.Series(p_val_dict)
outlier_df.name = "pval_outlier"

output_csv_file = os.path.join(config["output_file_name_list"]["OUTPUT_FOLDER"],
                               config["output_file_name_list"]["targeting_outlier_table"])

try:
    outlier_df.to_csv(output_csv_file, header=True)
    print(f"Successfully saved outlier p-values to '{output_csv_file}'")
except Exception as e:
    print(f"Error: Failed to save results to '{output_csv_file}'. Error: {e}")


### Filtering non-targeting gRNAs


non_target_gRNA_list = [item for item in list(gRNA_dict.keys())
                        if (item.startswith("non-targeting")) and
                        (len(gRNA_dict[item])>20)]
cell_id_nontarget_list = [gRNA_dict[key] for key in non_target_gRNA_list]

res = energy_distance_calc.pairwise_torch(pca_df,cell_id_nontarget_list,device,vardose=True)

# make Dataframe from results
pairwise_list = np.zeros((len(non_target_gRNA_list),
                          len(non_target_gRNA_list)
                         ))
for p1, p2, val in res:
    pairwise_list[p1,p2]=val
    pairwise_list[p2][p1]=val

df = pd.DataFrame(pairwise_list.copy(),
                  index=non_target_gRNA_list,
                  columns=non_target_gRNA_list)

df.index.name = "sgRNA"
df.columns.name = "sgRNA"
df.name = 'pairwise PCA distances'

sigmas = np.diag(df.values)
target_estats = 2 * df.values - sigmas - sigmas[:, np.newaxis]

kmeans = KMeans(n_clusters=3, random_state=0, n_init="auto").fit(target_estats)
label,count = np.unique(kmeans.labels_,return_counts=True)
print("k-means analysis of non-targeting")
print(label,count)

# Choose largest group as a "real" non-targeting background
largest_group_label = label[np.argmax(count)]
largest_group_ratio = np.round(np.sum(kmeans.labels_==largest_group_label)/len(kmeans.labels_),3)

print("largest group: ",largest_group_label)
print("largest group ratio: ",largest_group_ratio)


non_target_gRNA_name_df = pd.DataFrame(index=non_target_gRNA_list)
non_target_gRNA_name_df["pval_outlier"] = (kmeans.labels_==largest_group_label).astype(int)



#output the result
nt_output_csv_file = os.path.join(config["output_file_name_list"]["OUTPUT_FOLDER"],
                                  config["output_file_name_list"]["non_targeting_outlier_table"])

non_target_gRNA_name_df.to_csv(nt_output_csv_file)
