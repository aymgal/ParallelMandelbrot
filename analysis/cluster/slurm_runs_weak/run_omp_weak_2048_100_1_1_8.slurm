#!/bin/bash
#SBATCH --job-name=omp_weak
#SBATCH --account=phpc2017
#SBATCH --mail-user=aymeric.galan@epfl.ch --mail-type=FAIL

####SBATCH --ntasks=
####SBATCH --nodes=
####SBATCH --cpu_per_task=
####SBATCH --reservation=

module purge
module load gcc/5.3.0
module load mvapich2/2.2b

outfile=omp_weak_x_100_1_1_x.dat

export OMP_NUM_THREADS=8

for (( i = 1; i <= 10; i++ )) # for averaging compute time over 10 samples
do	
	srun ../../src/mandel_omp 2048 100 1 >> outputs/$outfile
done