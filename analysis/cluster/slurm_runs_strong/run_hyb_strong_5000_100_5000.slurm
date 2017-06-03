#!/bin/bash
#SBATCH --job-name=hyb_strong
#SBATCH --account=phpc2017
#SBATCH --mail-user=aymeric.galan@epfl.ch --mail-type=FAIL

####SBATCH --reservation=

module purge
module load gcc/5.3.0
module load mvapich2/2.2b

outfile=hyb_strong_5000_100_5000_x_4.dat

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

for (( i = 1; i <= 10; i++ )) # for averaging compute time over 10 samples
do	
	srun ../../src/mandel_hyb 5000 100 5000 >> outputs/$outfile
done