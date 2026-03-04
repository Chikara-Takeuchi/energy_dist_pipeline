#!/usr/bin/env python
# coding: utf-8

import sys
import os

import pandas as pd
import numpy as np
import scanpy as sc

from itertools import combinations
from sklearn.manifold import TSNE
from sklearn.cluster import AffinityPropagation

from tqdm import tqdm
import gc
import json
import torch

import util_functions
import energy_distance_calc

### Load PCA matric and gRNA matrix, and related information
json_fp = sys.argv[1]
json_fp_cluster = sys.argv[2]

with open(json_fp, 'r') as fp:
    config = json.load(fp)

with open(json_fp_cluster, 'r') as fp:
    config_clustering = json.load(fp)
    
print("--- Configuration Loaded ---")
print(json.dumps(config_clustering, indent=4))
print("--------------------------")

device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

(pca_df,gRNA_dict) = util_functions.load_files(config["input_data"]["h5ad_file"]["file_path"],
                                               config["input_data"]["sgRNA_file"]["file_path"],
                                               os.path.join(config["output_file_name_list"]["OUTPUT_FOLDER"],
                                                            config["output_file_name_list"]["pca_table"]),
                                               os.path.join(config["output_file_name_list"]["OUTPUT_FOLDER"],
                                                            config["output_file_name_list"]["gRNA_dict"]),
                                               obsm_key=config["input_data"]["h5ad_file"]["obsm_key"],
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

annotation_df = util_functions.load_annotation(config["input_data"]["annotation_file"]["file_path"])

gRNA_region_dict = util_functions.get_gRNA_region_dict(annotation_df,
                                                       gRNA_dict,
                                                       config["input_data"]["annotation_file"]["concatenate_key"])


#gRNA_region_clear_dict: basically gRNA_region_dict but outlier gRNAs are excluded
#cell_per_region_dict: cell IDs per each region
gRNA_region_clear_dict = {}

for key in gRNA_region_dict.keys():
    gRNA_list_tmp = [x for x in gRNA_region_dict[key] if x in clear_sgRNA_list]
    if len(gRNA_list_tmp)!=0:
        gRNA_region_clear_dict[key] = [x for x in gRNA_region_dict[key] if x in clear_sgRNA_list]

cell_per_region_dict = {}
for key in gRNA_region_clear_dict.keys():
    cell_list_tmp = [gRNA_dict[i] for i in gRNA_region_clear_dict[key]]
    cell_list_tmp = np.concatenate(cell_list_tmp)
    cell_per_region_dict[key] = np.unique(cell_list_tmp)



# Load pval_df from step2
pval_df = pd.read_csv(os.path.join(config["output_file_name_list"]["OUTPUT_FOLDER"],
                                   config["output_file_name_list"]["edist_pvalue_table"]),index_col=0)


### Calculate energy distances between perturbations
pval_df_sig = pval_df[(pval_df["pval_mean"]<config_clustering["cutoff"]["pval_cutoff"]) &
                      (pval_df["distance_mean"]>config_clustering["cutoff"]["distance_cutoff"]) &
                      (pval_df["type"]=="targeting")]

region_list_sig = np.unique(pval_df_sig.index)

cell_id_list_target = [cell_per_region_dict[key] for key in region_list_sig]


combi = list(combinations(range(len(cell_id_list_target)),2)) + \
        [(x,x) for x in range(len(cell_id_list_target))]


downsampling = config["aggregate"]["downsampling_maximum"]
res = []

# Main loop
pbar=tqdm(combi)
for target1_idx, target2_idx in pbar:
    pbar.set_postfix({"target1":region_list_sig[target1_idx],
                      "target2":region_list_sig[target2_idx],
                      "mode": "GPU"
                     })
    cell_test1 = cell_id_list_target[target1_idx]
    cell_test2 = cell_id_list_target[target2_idx]

    if len(cell_test1) > downsampling:
        cell_test1 = cell_test1[:downsampling]
    if len(cell_test2) > downsampling:
        cell_test2 = cell_test2[:downsampling]

    obs_edist = None
    mode = device # Assume GPU initially

    try:
        # Attempt GPU calculation
        obs_edist = energy_distance_calc.permutation_test(
            pca_df, cell_test1, cell_test2,
            device, 1, 1,return_permute=False
        )

    except Exception as e_gpu:
        print(f"GPU calculation failed: {e_gpu}. Attempting CPU fallback...")
        # Clean up GPU memory before CPU attempt
        gc.collect()
        if torch.cuda.is_available():
             torch.cuda.empty_cache()

        mode = "CPU" # Switch mode for postfix and calculation

        try:
            # Attempt CPU calculation
            obs_edist = energy_distance_calc.permutation_test(
                pca_df, cell_test1, cell_test2,
                "cpu", 1, 1,return_permute=False
            )

        except Exception as e_cpu:
            # Both GPU and CPU attempts failed
            print(f"CPU calculation also failed: {e_cpu}")
            print("Skipping energy distance calculation for this iteration.")
            # Return None to indicate failure


    res += [(region_list_sig[target1_idx],region_list_sig[target2_idx],obs_edist.item())]


#Convert dict format to dataframe
pairwise_dict=dict(zip(region_list_sig,[" "]*len(region_list_sig)))
for key in pairwise_dict.keys():
    pairwise_dict[key] = dict(zip(region_list_sig,[" "]*len(region_list_sig)))

for p1, p2, val in tqdm(res):
    pairwise_dict[p1][p2]=val
    pairwise_dict[p2][p1]=val
target_estats = pd.DataFrame(pairwise_dict,index=region_list_sig,columns=region_list_sig)

#Output to csv
target_estats.to_csv(os.path.join(config["output_file_name_list"]["OUTPUT_FOLDER"],
                                  config["output_file_name_list"]["edist_target_by_target_matrix"]))


### Clustering and embedding

#Embedding using tSNE
fit_method = TSNE(n_components=2, perplexity=2,n_iter=5000,random_state=1,
                  init="random",metric="precomputed")
embedding = fit_method.fit_transform(target_estats.copy())

total_edist_emb=pd.DataFrame(embedding,index=target_estats.index,columns=["x","y"]).reset_index()


# For now, only affinity propagation is supported
if config_clustering["clustering"]["method"]=="Affinity":
    clustering_method_emb = AffinityPropagation(random_state=0,
                                                convergence_iter=15,damping=0.50,
                                                affinity="precomputed"
                                               )
    cluster_info_emb = clustering_method_emb.fit(-target_estats)
else:
    clustering_method_emb = AffinityPropagation(random_state=0,
                                                convergence_iter=15,damping=0.50,
                                                affinity="precomputed"
                                               )
    cluster_info_emb = clustering_method_emb.fit(-target_estats)

total_edist_emb["cluster"] = cluster_info_emb.labels_

total_edist_emb.to_csv(os.path.join(config["output_file_name_list"]["OUTPUT_FOLDER"],
                                  config["output_file_name_list"]["edist_embedding_info"]))
