#!/usr/bin/bash
#
#SBATCH -J edist_pipeline_step0       # Job name
#SBATCH -N 1                          # Total number of nodes requested (16 cores/node)
#SBATCH -t 24:00:00                   # Run time (hh:mm:ss) - 20 hrs limit
#SBATCH -p 256GBv2
#SBATCH -o run_output_step0.out
#SBATCH -e run_output_step0.err


#Note: Use singularity 4.1.0 for this code. for UTSW environment, run "module load singularity/4.1.0" to activate singularity. 
# module load singularity/4.1.0
#Define the path to the config file and bin directory
CONTAINER_PATH="/project/GCRB/Hon_lab/s223695/Data_project/Perturb_seq_edist_pipeline/container/edist_pipeline_v01.sif"

singularity exec --nv ${CONTAINER_PATH} bash ./run_step0_base.sh