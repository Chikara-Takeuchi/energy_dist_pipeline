#!/usr/bin/bash
#
#SBATCH -Jedist_pipeline_step1_2      # Job name
#SBATCH -N 1                          # Total number of nodes requested (16 cores/node)
#SBATCH -t 24:00:00                   # Run time (hh:mm:ss) - 20 hrs limit
#SBATCH -p GPUv100s
#SBATCH -o run_output_step12.out
#SBATCH -e run_output_step12.err
#SBATCH --gres=gpu:1

#Note: Use singularity 4.1.0 for this code. for UTSW environment, run "module load singularity/4.1.0" to activate singularity. 
module load singularity/4.1.0

#Define the path to the config file and bin directory
CONTAINER_PATH="/project/GCRB/Hon_lab/s223695/Data_project/Perturb_seq_edist_pipeline/container/edist_pipeline_v01.sif"
CONFIG_PATH="/project/GCRB/Hon_lab/s223695/Data_project/Perturb_seq_edist_pipeline/config.json"
BIN_PATH="/project/GCRB/Hon_lab/s223695/Data_project/Perturb_seq_edist_pipeline/bin"

#Note: When you run this code, please make sure that /pipeline_output/annotation_file_table.csv (or a defined name in config file) exists.

echo "[Step1] Filtering outlier gRNAs"
singularity exec --nv ${CONTAINER_PATH} python ${BIN_PATH}/1_filtereing_gRNA.py ${CONFIG_PATH}

echo "[Step2] calculate energy distance between targets and non-targeting"
singularity exec --nv ${CONTAINER_PATH} python ${BIN_PATH}/2_e_distance_nontargeting.py ${CONFIG_PATH}

echo "[Step2_1] visualize results of energy distance analysis"
singularity exec --nv ${CONTAINER_PATH} python ${BIN_PATH}/2_1_Plot_figure.py ${CONFIG_PATH}

echo "All steps completed successfully."