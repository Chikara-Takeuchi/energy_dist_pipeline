#!/usr/bin/bash
#
#SBATCH -J edist_pipeline_DRYRUN      # Job name
#SBATCH -N 1                          # Total number of nodes requested (16 cores/node)
#SBATCH -t 24:00:00                   # Run time (hh:mm:ss) - 20 hrs limit
#SBATCH -p GPUp4
#SBATCH -o run_output_dryrun.out
#SBATCH -e run_output_dryrun.err
#SBATCH --gres=gpu:1


#Define the path to the config file and bin directory
CONTAINER_PATH="/project/GCRB/Hon_lab/s223695/Data_project/Perturb_seq_edist_pipeline/container/edist_pipeline_v01.sif"

#Note: When you run this code, please make sure that /pipeline_output/annotation_file_table.csv (or a defined name in config file) exists.
singularity exec --nv ${CONTAINER_PATH} bash ./DRYRUN_FILES/run_step1_2_base_DRYRUN.sh
singularity exec --nv ${CONTAINER_PATH} bash ./DRYRUN_FILES/run_step3_base_DRYRUN.sh