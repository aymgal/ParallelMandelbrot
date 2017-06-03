#!/bin/bash

#----------------------------- strong scaling --------------------------------#

sdir=slurm_runs_strong

N=5000
n_iter=100
n_row=$N
n_thread=1

outfile=mpi_strong_${N}_${n_iter}_${n_row}_x_${n_thread}.dat
echo "Output file name : " $outfile

mkdir -p outputs

rm -f outputs/$outfile
rm -f out_*.pgm
rm -f out_*.bmp

cd ../../src 			# go there to compile
make mpi
cd ../analysis/cluster	# come back here to run

for n_proc in 2 4 8 16 32 64
do
	echo "Going to run several times mandel_mpi, with parameters N = $N, \
max iter = $n_iter, rows = $n_row, procs = $n_proc, threads = $n_thread)"

	sbatch -n ${n_proc} ${sdir}/run_mpi_strong_${N}_${n_iter}_${n_row}_${n_thread}.slurm
done

#-----------------------------------------------------------------------------#
