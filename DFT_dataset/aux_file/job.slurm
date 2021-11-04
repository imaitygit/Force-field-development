#!/bin/bash --login

#SBATCH --job-name=kc
#SBATCH --nodes=1
#SBATCH --tasks-per-node=48
#SBATCH --cpus-per-task=1
#SBATCH --time=00:20:00

# Replace [budget code] below with your project code (e.g. t01)
#SBATCH --account=e05-nano-lis
#SBATCH --partition=standard
#SBATCH --qos=short
#SBATCH --reservation=shortqos


# Ensure that the environment is set correctly for the job
module load epcc-job-env

# Load quantum espresso
module load quantum-espresso/6.6

# Load the xthi module to get access to the xthi program
#module load xthi

# srun launches the parallel program based on the SBATCH options
# The SBATCH options are provided above
# This is an "only MPI" job.
srun --distribution=block:block --hint=nomultithread --unbuffered pw.x -in pw.in >& pw.out
