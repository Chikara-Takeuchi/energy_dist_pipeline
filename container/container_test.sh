#!/bin/bash

# Name for the job that will be visible in the job queue and accounting tools.
#SBATCH --job-name NewJob
#SBATCH -p GPU       # partition (queue)
#SBATCH --gres=gpu:1
# Number of nodes required to run this job
#SBATCH -N 1
#SBATCH -t 0-2:0:0

# The standard output and errors from commands will be written to these files.
# %j in the filename will be replace with the job number when it is submitted.
#SBATCH -o job_%j.out
#SBATCH -e job_%j.err

# Send an email when the job status changes, to the specfied address.
#SBATCH --mail-type ALL
#SBATCH --mail-user Chikara.Takeuchi@utsouthwestern.edu

module load singularity/3.9.9
singularity exec --nv ./edist_pipeline_v01.sif nvidia-smi



# END OF SCRIPT
