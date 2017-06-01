#!/bin/bash
#SBATCH --ntasks-per-node=1
#SBATCH --job-name=hyb_strong
#SBATCH --account=phpc2017
#SBATCH --mail-user=aymeric.galan@epfl.ch --mail-type=FAIL

####SBATCH --ntasks=
####SBATCH --nodes=
####SBATCH --cpu_per_task=
####SBATCH --reservation=

module purge
module load gcc/5.3.0
module load mvapich2/2.2b

outfile=hyb_strong_10000_100_1_2_x.dat

export OMP_NUM_THREAD=1

for (( i = 1; i <= 10; i++ )) # for averaging compute time over 10 samples
do	
	srun ../../src/mandel_hyb 10000 100 1 >> outputs/$outfile
done