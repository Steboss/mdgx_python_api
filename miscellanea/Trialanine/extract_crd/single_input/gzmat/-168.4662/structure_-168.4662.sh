#!/bin/bash
#SBATCH -o output-%A-%a.out
#SBATCH -p serial -n 8
#SBATCH --time 48:00:00

source /etc/profile.d/module.sh
export OMP_NUM_THREADS=1
g09 < structure_-168.4662.gzmat > structure_-168.4662.gout

wait
