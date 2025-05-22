#!/usr/bin/bash

#Define the path to the config file and bin directory
CONFIG_PATH="/project/GCRB/Hon_lab/s223695/Data_project/Perturb_seq_edist_pipeline/config.json"
BIN_PATH="/project/GCRB/Hon_lab/s223695/Data_project/Perturb_seq_edist_pipeline/bin"

# Load the conda module
echo "[Step0] preprocessing annotation files"
python ${BIN_PATH}/0_preprocess.py ${CONFIG_PATH}

echo "All steps completed successfully."