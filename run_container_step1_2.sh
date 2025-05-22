#!/usr/bin/bash
#
#SBATCH -Jedist_pipeline_step1_2      # Job name
#SBATCH -N 1                          # Total number of nodes requested (16 cores/node)
#SBATCH -t 24:00:00                   # Run time (hh:mm:ss) - 20 hrs limit
#SBATCH -p GPUv100s
#SBATCH -o run_output_step12.out
#SBATCH -e run_output_step12.err
#SBATCH --gres=gpu:1


#Define the path to the config file and bin directory
CONTAINER_PATH="/project/GCRB/Hon_lab/s223695/Data_project/Perturb_seq_edist_pipeline/container/edist_pipeline_v01.sif"

#Note: When you run this code, please make sure that /pipeline_output/annotation_file_table.csv (or a defined name in config file) exists.
singularity exec --nv ${CONTAINER_PATH} bash ./run_step1_2_base.sh