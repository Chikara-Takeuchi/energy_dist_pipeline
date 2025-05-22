#!/usr/bin/bash
#
#SBATCH -J edist_pipeline_step3       # Job name
#SBATCH -N 1                          # Total number of nodes requested (16 cores/node)
#SBATCH -t 24:00:00                   # Run time (hh:mm:ss) - 20 hrs limit
#SBATCH -p GPUv100s
#SBATCH -o run_output_step3.out
#SBATCH -e run_output_step3.err
#SBATCH --gres=gpu:1


#Define the path to the config file and bin directory
CONTAINER_PATH="/project/GCRB/Hon_lab/s223695/Data_project/Perturb_seq_edist_pipeline/container/edist_pipeline_v01.sif"

singularity exec --nv ${CONTAINER_PATH} bash ./run_step3_base.sh