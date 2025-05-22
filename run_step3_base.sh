#!/usr/bin/bash

#Define the path to the config file and bin directory
CONFIG_PATH="/project/GCRB/Hon_lab/s223695/Data_project/Perturb_seq_edist_pipeline/config.json"
CONFIG_CLUSTER_PATH="/project/GCRB/Hon_lab/s223695/Data_project/Perturb_seq_edist_pipeline/config_clustering.json"

BIN_PATH="/project/GCRB/Hon_lab/s223695/Data_project/Perturb_seq_edist_pipeline/bin"

echo "[Step3] Calculate association between perturbations"
python ${BIN_PATH}/3_e_distance_among_regions.py \
  ${CONFIG_PATH} ${CONFIG_CLUSTER_PATH}

echo "All steps completed successfully."
