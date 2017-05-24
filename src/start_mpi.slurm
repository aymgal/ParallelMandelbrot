#!/bin/bash
#SBATCH --reservation phpc2017
#SBATCH --account phpc2017
#SBATCH --ntasks 1

mpirun ./mandel_mpi 1000 1000