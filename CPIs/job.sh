#!/bin/bash 
#SBATCH -N 1
#SBATCH -p normal
#SBATCH --gres=gpu:2
#SBATCH -n 6
#SBATCH -o meld.log

# Other parameters for slurm may be added 
# For MELD, the number of tasks much match the number of allocated GPUs!

python CPI_Adapt.py

mpirun -np 2 launch_remd