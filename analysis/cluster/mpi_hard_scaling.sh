#!/bin/bash

#------------------------------ hard scaling ---------------------------------#

N=1000
n_iter=100
n_row=100
n_threads=0 	# because unused in the case of non-OpenMP code

outfile=mpi_hard_${N}_${n_iter}_${n_row}_x_${n_threads}.dat
echo "Output file name : " $outfile

rm -f outputs/$outfile
rm -f out_*.pgm
rm -f out_*.bmp

cd ../../src 			# go there to compile
make mpi
cd ../analysis/slurm	# come back here to run

n_proc_min=2
n_proc_max=32
n_proc=$n_proc_min

while [[ $n_proc -le $n_proc_max ]]
do
	echo "Running several times mandel_mpi, with parameters N = $N, \
max iter = $n_iter, rows = $n_row, procs = $n_proc, threads = $n_threads)"

	sbatch -n ${n_proc} run_mpi_hard.slurm

	((n_proc = n_proc * 2)) # grows number of processes, at fixed system size
done

#-----------------------------------------------------------------------------#
