#!/bin/bash
#SBATCH --job-name=mpi_proc
####SBATCH --nodes=32
#SBATCH --ntasks=32
####SBATCH --cpu_per_task=1
#SBATCH --account=phpc2017
####SBATCH --reservation=phpc2017
#SBATCH --mail-user=aymeric.galan@epfl.ch --mail-type=ALL

rm -f outputs/mpi_var_proc.dat
rm -f out_*.pgm
rm -f out_*.bmp

cd ../../src 			# go there to compile
make mpi
cd ../analysis/slurm	# come back here to run

#------------------------------ hard scaling ---------------------------------#

N=1000
n_iter=100
n_row=100
n_threads=0 	# because unused in the case of non-OpenMP code

outfile=mpi_hard_${N}_${n_iter}_${n_row}_x_${n_threads}.dat
echo "Output file name : " $outfile

n_proc_min=2
n_proc_max=32
n_proc=$n_proc_min

while [[ $n_proc -le $n_proc_max ]]
do
	for (( i = 1; i <= 10; i++ )) # for averaging compute time over 10 samples
	do
		echo "Running mandel_mpi (N = $N, max iter = $n_iter, \
rows = $n_row, procs = $n_proc, threads = $n_threads)"

		mpirun -np $n_proc ../../src/mandel_mpi $N $n_iter $n_row \
>> outputs/$outfile
	done
	((n_proc = n_proc * 2))
done
#-----------------------------------------------------------------------------#
