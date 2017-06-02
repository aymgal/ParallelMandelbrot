#!/bin/bash
#SBATCH --job-name=hyb_strong_thr
#SBATCH --ntasks=2		# number of MPI ranks
#SBATCH --cpus-per-task=2	# number of threads/rank
#SBATCH --ntasks-per-node=1
#SBATCH --ntasks-per-socket=1
#SBATCH --account=phpc2017
#SBATCH --mail-user=aymeric.galan@epfl.ch --mail-type=FAIL

####SBATCH --reservation=

module purge
module load gcc/5.3.0
module load mvapich2/2.2b

outfile=hyb_strong_1000_100_1_2_x.dat

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

for (( i = 1; i <= 10; i++ )) # for averaging compute time over 10 samples
do	
	srun ../../src/mandel_hyb 1000 100 1 >> outputs/$outfile
done