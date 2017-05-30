#!/bin/bash
#SBATCH --job-name=mpi_soft
#SBATCH --account=phpc2017
#SBATCH --mail-user=aymeric.galan@epfl.ch --mail-type=FAIL

####SBATCH --ntasks=
####SBATCH --nodes=
####SBATCH --cpu_per_task=
####SBATCH --reservation=

outfile=mpi_weak_x_${n_iter}_${n_row}_x_${n_threads}.dat

for (( i = 1; i <= 10; i++ )) # for averaging compute time over 10 samples
do	
	srun ../../src/mandel_mpi 512 $n_iter $n_row >> outputs/$outfile
done