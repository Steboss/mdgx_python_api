#!/bin/bash
#SBATCH -o output-%A-%a.out
#SBATCH -p serial -n 8
#SBATCH --time 48:00:00

source /etc/profile.d/module.sh
export OMP_NUM_THREADS=1
g09 < structure_340.0.gzmat > structure_340.0.gout

wait
