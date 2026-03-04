#!/usr/bin/bash
#
#SBATCH -J edist_pipeline_step1_2      # Job name
#SBATCH -N 1                          # Total number of nodes requested (16 cores/node)
#SBATCH -t 24:00:00                   # Run time (hh:mm:ss) - 20 hrs limit
#SBATCH -p GPUv100s
#SBATCH -o run_output_step12.out
#SBATCH -e run_output_step12.err
#SBATCH --gres=gpu:1

module load apptainer
set -euo pipefail

nvidia-smi

#Define the path to the config file and bin directory
CONTAINER_PATH="/project/GCRB/Hon_lab/s223695/Data_project/Perturb_seq_edist_pipeline/container/edist_pipeline_v2.sif"
CONFIG_PATH="/project/GCRB/Hon_lab/s223695/Data_project/Perturb_seq_edist_pipeline/stable_version/config.json"
BIN_PATH="/project/GCRB/Hon_lab/s223695/Data_project/Perturb_seq_edist_pipeline/stable_version/bin"

#Note: When you run this code, please make sure that /pipeline_output/annotation_file_table.csv (or a defined name in config file) exists.

echo "[Step1] Filtering outlier gRNAs"
apptainer exec -B /project:/project --nv ${CONTAINER_PATH} python ${BIN_PATH}/1_filtering_gRNA.py ${CONFIG_PATH}

echo "[Step2] calculate energy distance between targets and non-targeting"
apptainer exec -B /project:/project --nv ${CONTAINER_PATH} python ${BIN_PATH}/2_e_distance_nontargeting.py ${CONFIG_PATH}

echo "[Step2_1] visualize results of energy distance analysis"
apptainer exec -B /project:/project --nv ${CONTAINER_PATH} python ${BIN_PATH}/2_1_Plot_figure.py ${CONFIG_PATH}

echo "All steps completed successfully."
