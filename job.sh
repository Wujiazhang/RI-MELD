#!/bin/bash 
#SBATCH -N 1
#SBATCH -p normal
#SBATCH --gres=gpu:3
#SBATCH -n 6
#SBATCH -o meld.log

# Other parameters for slurm may be added 
# For MELD, the number of tasks much match the number of allocated GPUs!

python CPI_adapt.py

mpirun -np 3 launch_remd