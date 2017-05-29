#!/bin/bash
#SBATCH --job-name=mpi_N
####SBATCH --nodes=32
#SBATCH --ntasks=8
####SBATCH --cpu_per_task=1
#SBATCH --account=phpc2017
####SBATCH --reservation=phpc2017
#SBATCH --mail-user=aymeric.galan@epfl.ch --mail-type=ALL

rm -f outputs/mpi_var_N.dat
rm -f out_*.pgm
rm -f out_*.bmp

cd ../../src 			# go there to compile
make mpi
cd ../analysis/slurm	# come back here to run

#------------------------- variable : problem size ---------------------------#

n_iter=100
n_row=100
n_proc=8
n_threads=0	# because unused in the case of non-OpenMP code

outfile=mpi_N_x_${n_iter}_${n_row}_${n_proc}_${n_threads}.dat
echo "Output file name : " $outfile

N_min=64
N_max=8192
N=$N_min

while [[ $N -le $N_max ]]
do
	for (( i = 1; i <= 10; i++ )) # for averaging compute time over 10 samples
	do
		echo "Running mandel_mpi (N = $N, max iter = $n_iter, \
rows = $n_row, procs = $n_proc, threads = $n_threads)"
		mpirun -np $n_proc ../../src/mandel_mpi $N $n_iter $n_row \
>> outputs/$outfile
	done
	((N = N * 2))
done

#-----------------------------------------------------------------------------#

