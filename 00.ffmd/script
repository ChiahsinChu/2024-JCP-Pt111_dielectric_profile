#!/bin/bash

#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=64
#SBATCH --partition=small_s
#SBATCH --job-name=cp2k
#SBATCH --output=%j.stdout 
#SBATCH --error=%j.stderr 
#SBATCH --mem=4G

##############################################
#          Software Envrironment             #
##############################################
module purge

module load gcc/9.3 intel/2020.2
module load cp2k/2022.1-intel-2020
##############################################
#               Run job                      #
##############################################
export OMP_NUM_THREADS=1
mpirun cp2k.psmp -i input.inp > output
