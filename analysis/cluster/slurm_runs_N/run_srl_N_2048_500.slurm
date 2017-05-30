#!/bin/bash
#SBATCH --job-name=srl_N
#SBATCH --account=phpc2017
#SBATCH --mail-user=aymeric.galan@epfl.ch --mail-type=FAIL

####SBATCH --ntasks=
####SBATCH --nodes=
####SBATCH --cpu_per_task=
####SBATCH --reservation=

module purge
module load gcc/5.3.0
module load mvapich2/2.2b

outfile=srl_N_x_100_1_1_1.dat

for (( i = 1; i <= 10; i++ )) # for averaging compute time over 10 samples
do	
	srun ../../src/mandel_srl 2048 500 1 >> outputs/$outfile
done