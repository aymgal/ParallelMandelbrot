#!/bin/bash
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus_per_task 1
#SBATCH --reservation phpc2017
#SBATCH --account phpc2017

mpirun ./mandel_mpi 1000 1000