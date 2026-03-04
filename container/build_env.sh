#!/usr/bin/bash
#
#SBATCH -J build_env      # Job name
#SBATCH -N 1                          # Total number of nodes requested (16 cores/node)
#SBATCH -t 24:00:00                   # Run time (hh:mm:ss) - 20 hrs limit
#SBATCH -p 256GBv2
#SBATCH -o ./build_log.out
#SBATCH -e ./build_log.err
uv lock
podman build -t e-dist_pipeline .
