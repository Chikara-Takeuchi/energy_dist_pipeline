#!/usr/bin/env python
# coding: utf-8

# Note: please modify the file path so that script runs in your environment
#
#
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os
import json

import util_functions
from tqdm import tqdm

# Config start
gRNA_ref_file_path="/project/GCRB/Hon_lab/s223695/Data_project/Perturb_seq_edist_pipeline/ref/Hon_sgRNA_index_dacc_annot_reference.csv"
neg_control_file_path="/project/GCRB/Hon_lab/s223695/Data_project/Perturb_seq_edist_pipeline/ref/negative_controls.tsv"
non_target_file_path="/project/GCRB/Hon_lab/s223695/Data_project/Perturb_seq_edist_pipeline/ref/non_targeting.tsv"

pos_control_gene_list = ["CD81","CD151","CD55","CD29","B2M","AARS","POLR1D","DNAJC19","MALAT1","NGFRP1","TFRC"]

json_fp = "./config.json"

# Config end
with open(json_fp, 'r') as fp:
    config = json.load(fp)

output_folder = config["output_file_name_list"]["OUTPUT_FOLDER"]

if os.path.exists(output_folder)==False:
    print("generate folder for figure:",output_folder)
    os.mkdir(output_folder)
else:
    print("Folder already exist",output_folder)
    
gRNA_ref_df = pd.read_csv(gRNA_ref_file_path,sep="\t")

neg_control_df = \
    pd.read_csv(neg_control_file_path,sep="\t",index_col=0)
non_target_df = \
    pd.read_csv(non_target_file_path,sep="\t",index_col=0)

neg_control_name = neg_control_df.index.tolist()
pos_control_name = pos_control_gene_list
non_target_name = non_target_df.index.tolist()

def detect_source(target_gRNA_name):
    target_gene = util_functions.extract_gene_name(target_gRNA_name)
    if target_gene in neg_control_name:
        return "negative control"
    elif target_gene in pos_control_name:
        return "positive control"
    elif target_gene=="non-targeting":
        return "non-targeting"
    else:
        return "targeting"


gRNA_ref_df["target_transcript_name"] = gRNA_ref_df["protospacer_ID"].apply(util_functions.extract_transcript_name)
gRNA_ref_df["source"] = gRNA_ref_df["protospacer_ID"].apply(detect_source)
gRNA_ref_df["target_gene_name"] = gRNA_ref_df["intended_target_name"].copy()

print(np.unique(gRNA_ref_df["source"],return_counts=True))


gRNA_ref_df_output = gRNA_ref_df.loc[:,["protospacer_ID","target_transcript_name","source","protospacer"]]

gRNA_ref_df_output.columns = ["guide_id","intended_target_name","type","spacer"]

gRNA_ref_df_output.to_csv(os.path.join(config["output_file_name_list"]["OUTPUT_FOLDER"],
                                       config["output_file_name_list"]["annotation_file"]))



