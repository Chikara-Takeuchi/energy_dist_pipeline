#!/usr/bin/bash

#Define the path to the config file and bin directory
CONFIG_PATH="/project/GCRB/Hon_lab/s223695/Data_project/Perturb_seq_edist_pipeline/config.json"
BIN_PATH="/project/GCRB/Hon_lab/s223695/Data_project/Perturb_seq_edist_pipeline/bin"

# Load the conda module

echo "[Step1] Filtering outlier gRNAs"
python ${BIN_PATH}/1_filtereing_gRNA.py ${CONFIG_PATH}

echo "[Step2] calculate energy distance between targets and non-targeting"
python ${BIN_PATH}/2_e_distance_nontargeting.py ${CONFIG_PATH}

echo "[Step2_1] visualize results of energy distance analysis"
python ${BIN_PATH}/2_1_Plot_figure.py ${CONFIG_PATH}

echo "All steps completed successfully."