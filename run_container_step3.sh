#!/usr/bin/bash
#
#SBATCH -J edist_pipeline_step3       # Job name
#SBATCH -N 1                          # Total number of nodes requested (16 cores/node)
#SBATCH -t 24:00:00                   # Run time (hh:mm:ss) - 20 hrs limit
#SBATCH -p GPUv100s
#SBATCH -o run_output_step3.out
#SBATCH -e run_output_step3.err
#SBATCH --gres=gpu:1

#Note: Use singularity 4.1.0 for this code. for UTSW environment, run "module load singularity/4.1.0" to activate singularity. 
module load singularity/4.1.0

#Define the path to the config file and bin directory
CONFIG_PATH="/project/GCRB/Hon_lab/s223695/Data_project/Perturb_seq_edist_pipeline/config.json"
CONFIG_CLUSTER_PATH="/project/GCRB/Hon_lab/s223695/Data_project/Perturb_seq_edist_pipeline/config_clustering.json"
BIN_PATH="/project/GCRB/Hon_lab/s223695/Data_project/Perturb_seq_edist_pipeline/bin"

# Load the conda module
source activate scanpy_gpu
echo "[Step3] Calculate association between perturbations"
singularity exec --nv ${CONTAINER_PATH}　python ${BIN_PATH}/3_e_distance_among_regions.py \
  ${CONFIG_PATH} ${CONFIG_CLUSTER_PATH}

echo "All steps completed successfully."