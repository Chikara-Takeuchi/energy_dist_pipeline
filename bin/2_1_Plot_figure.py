#!/usr/bin/env python
# coding: utf-8

import sys
import os

import pandas as pd
import numpy as np

import seaborn as sns
import matplotlib.pyplot as plt
plt.rcParams.update({'axes.labelsize' : 'large',
                     'pdf.fonttype' : 42
                    })

from tqdm import tqdm

import gc
import warnings
import time
import pickle
import json

import util_functions


### Load congigurations and related files
json_fp = sys.argv[1]
with open(json_fp, 'r') as fp:
    config = json.load(fp)


figure_folder = os.path.join(config["output_file_name_list"]["OUTPUT_FOLDER"],
                             "image")
if os.path.exists(figure_folder)==False:
    print("generate folder for figure:",figure_folder)
    os.mkdir(figure_folder)
else:
    print("Folder already exist",figure_folder)


sgRNA_outlier_df = pd.read_csv(os.path.join(config["output_file_name_list"]["OUTPUT_FOLDER"],
                                                  config["output_file_name_list"]["targeting_outlier_table"]),
                                     index_col=0)

nontargeting_outlier_df = pd.read_csv(os.path.join(config["output_file_name_list"]["OUTPUT_FOLDER"],
                                                         config["output_file_name_list"]["non_targeting_outlier_table"]),
                                            index_col=0)


annotation_df = util_functions.load_annotation(config["input_data"]["annotation_file"]["file_path"])


### plot number of (1) gRNAs removed in the outlier analysis and (2) gRNAs remained

grouping_name = config["input_data"]["annotation_file"]["concatenate_key"]

target_transcript_name_list = np.unique(annotation_df[annotation_df["type"]!="non-targeting"][grouping_name])
sgRNA_outlier_dict = sgRNA_outlier_df.to_dict()["pval_outlier"]
gRNA_stat_dict = {}

for target_name in tqdm(target_transcript_name_list):
    gRNA_names_all = annotation_df[annotation_df[grouping_name]==target_name]["guide_id"].values
    gRNA_stat_dict[target_name] = {}

    gRNA_dropped = 0
    gRNA_outlier = 0
    gRNA_remained = 0

    for gRNA_name in gRNA_names_all:
        if gRNA_name in sgRNA_outlier_dict.keys():
            if sgRNA_outlier_dict[gRNA_name]<0.05:
                gRNA_outlier = gRNA_outlier+1
            else:
                gRNA_remained = gRNA_remained+1
        else:
            gRNA_dropped = gRNA_dropped+1

    gRNA_stat_dict[target_name]["dropped"] = gRNA_dropped
    gRNA_stat_dict[target_name]["outlier"] = gRNA_outlier
    gRNA_stat_dict[target_name]["remained"] = gRNA_remained

gRNA_stat_df = pd.DataFrame(gRNA_stat_dict).T

fig,ax = plt.subplots(3,1,figsize=(10,6))
sns.histplot(data=gRNA_stat_df,x="dropped",ax=ax[0])
sns.histplot(data=gRNA_stat_df,x="outlier",ax=ax[1])
sns.histplot(data=gRNA_stat_df,x="remained",ax=ax[2])

plt.tight_layout()
plt.savefig(os.path.join(figure_folder,"gRNA_stat.pdf"))


#### Plot distirbution of energy distance to find cutoff

pval_edit_df = pd.read_csv(os.path.join(config["output_file_name_list"]["OUTPUT_FOLDER"],
                                        config["output_file_name_list"]["edist_pvalue_table"]),
                           index_col=0)
fig,ax = plt.subplots(1,1,figsize=(8,6))
sns.scatterplot(data=pval_edit_df,x="distance_mean_log",y="pval_mean_log")
sns.kdeplot(data=pval_edit_df,x="distance_mean_log",y="pval_mean_log",color="black",alpha=0.4)
plt.savefig(os.path.join(figure_folder,"e-dist_distribution.pdf"))


#### Find best cutoff
max_pval_round = int(np.round(np.max(pval_edit_df["pval_mean_log"])))

pval_list = [1/np.power(10,x) for x in list(range(max_pval_round+1))]
dist_val_list = [0,5,10,50,100]

# Count number of regions in each type
type_list,type_count = np.unique(pval_edit_df["type"],return_counts=True)
print("Number of regions in each type:")
for t,c in zip(type_list,type_count):
    print (t,":",c) 
if  "negative control" not in type_list:
    print("Warning: 'negative control' not found in the type list. Please check your data.")

cutoff_df = pd.DataFrame(index=pval_list,columns=dist_val_list)

for pval_cutoff in pval_list:
    for dist_cutoff in dist_val_list:
        selected_df = pval_edit_df[
                        (pval_edit_df["distance_mean"]>dist_cutoff) &
                        (pval_edit_df["pval_mean"]<pval_cutoff)]
        cutoff_df.loc[pval_cutoff,dist_cutoff] = selected_df.shape[0]

cutoff_df_neg = pd.DataFrame(index=pval_list,columns=dist_val_list)

for pval_cutoff in pval_list:
    for dist_cutoff in dist_val_list:
        selected_df = pval_edit_df[
                        (pval_edit_df["type"]=="negative control") &
                        (pval_edit_df["distance_mean"]>dist_cutoff) &
                        (pval_edit_df["pval_mean"]<pval_cutoff)]
        cutoff_df_neg.loc[pval_cutoff,dist_cutoff] = selected_df.shape[0]


# Find cutoff for all regions (targets, neg_control, pos_control)
fig,ax = plt.subplots(figsize=(6,5))
sns.heatmap(cutoff_df.astype(float),annot=True,fmt='.0f',ax=ax)
ax.set_xlabel("distance")
ax.set_ylabel("p-value")
plt.savefig(os.path.join(figure_folder,"e-dist_cutoff_value.pdf"))


# Find cutoff for all negative controls (neg_control)
fig,ax = plt.subplots(figsize=(6,5))
sns.heatmap(cutoff_df_neg.astype(float),annot=True,fmt='.0f',ax=ax)
ax.set_xlabel("distance")
ax.set_ylabel("p-value")
plt.savefig(os.path.join(figure_folder,"e-dist_cutoff_value_NEG_CONTROL.pdf"))
