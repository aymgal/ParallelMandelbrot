#!/bin/bash
#SBATCH --job-name=mpi_strong
#SBATCH --account=phpc2017
#SBATCH --mail-user=aymeric.galan@epfl.ch --mail-type=FAIL

####SBATCH --ntasks=
####SBATCH --nodes=
####SBATCH --cpu_per_task=
####SBATCH --reservation=

module purge
module load gcc/5.3.0
module load mvapich2/2.2b

N=10000
n_iter=100
n_row=100
n_threads=0 	# because unused in the case of non-OpenMP code

outfile=mpi_strong_${N}_${n_iter}_${n_row}_x_${n_threads}.dat

for (( i = 1; i <= 10; i++ )) # for averaging compute time over 10 samples
do	
	srun ../../src/mandel_mpi $N $n_iter $n_row >> ../outputs/$outfile
done