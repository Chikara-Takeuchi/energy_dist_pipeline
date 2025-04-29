#!/usr/bin/bash
## srun -N 1 -n 1 -p super --pty /bin/bash

module load singularity/3.9.9
echo "[login with podman]"
podman login git.biohpc.swmed.edu:5050

echo "[login with singularity]"
singularity remote login --username s223695 docker://git.biohpc.swmed.edu:5050

podman build -t git.biohpc.swmed.edu:5050/s223695/container_loc/edist_pipeline:v01 --format docker .

python3 podman.py push git.biohpc.swmed.edu:5050/s223695/container_loc/edist_pipeline:v01

singularity pull docker://git.biohpc.swmed.edu:5050/s223695/container_loc/edist_pipeline:v01
