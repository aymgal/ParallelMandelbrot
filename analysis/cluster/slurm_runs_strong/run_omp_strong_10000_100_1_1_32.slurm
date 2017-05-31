#!/bin/bash
#SBATCH --job-name=omp_strong
#SBATCH --account=phpc2017
#SBATCH --mail-user=aymeric.galan@epfl.ch --mail-type=FAIL

####SBATCH --ntasks=
####SBATCH --nodes=
####SBATCH --cpu_per_task=
####SBATCH --reservation=

module purge
module load gcc/5.3.0

outfile=omp_strong_10000_100_1_1_x.dat

export OMP_NUM_THREADS=32

for (( i = 1; i <= 10; i++ )) # for averaging compute time over 10 samples
do	
	srun ../../src/mandel_omp 10000 100 1 >> outputs/$outfile
done