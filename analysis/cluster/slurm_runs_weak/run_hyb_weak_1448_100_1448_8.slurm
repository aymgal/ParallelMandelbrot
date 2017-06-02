#!/bin/bash
#SBATCH --job-name=hyb_weak
#SBATCH --account=phpc2017
#SBATCH --mail-user=aymeric.galan@epfl.ch --mail-type=FAIL

####SBATCH --ntasks=
####SBATCH --nodes=
####SBATCH --cpu_per_task=
####SBATCH --ntasks-per-node=
####SBATCH --ntasks-per-socket=
####SBATCH --reservation=

module purge
module load gcc/5.3.0
module load mvapich2/2.2b

outfile=hyb_weak_x_100_x_x_8.dat

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

for (( i = 1; i <= 10; i++ )) # for averaging compute time over 10 samples
do	
	srun ../../src/mandel_hyb 1448 100 1448 >> outputs/$outfile
done