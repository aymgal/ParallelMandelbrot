#!/bin/bash
#SBATCH --job-name=hyb_N
#SBATCH --account=phpc2017
#SBATCH --mail-user=aymeric.galan@epfl.ch --mail-type=FAIL

####SBATCH --ntasks=
####SBATCH --nodes=
####SBATCH --cpu_per_task=
####SBATCH --reservation=

module purge
module load gcc/5.3.0
module load mvapich2/2.2b

outfile=hyb_N_x_100_x_16_4.dat

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

for (( i = 1; i <= 10; i++ )) # for averaging compute time over 10 samples
do	
	srun ../../src/mandel_hyb 4096 100 4096 >> outputs/$outfile
done