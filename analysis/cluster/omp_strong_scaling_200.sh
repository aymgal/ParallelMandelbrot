#!/bin/bash

#----------------------------- strong scaling --------------------------------#

sdir=slurm_runs_strong

N=10000
n_iter=200
n_row=1
n_proc=1

outfile=omp_strong_${N}_${n_iter}_${n_row}_${n_proc}_x.dat
echo "Output file name : " $outfile

mkdir -p outputs

rm -f outputs/$outfile
rm -f out_*.pgm
rm -f out_*.bmp

cd ../../src 			# go there to compile
make omp
cd ../analysis/cluster	# come back here to run

for n_thread in 1 2 4 8 16 32 64
do
	echo "Going to run several times mandel_omp, with parameters N = $N, \
max iter = $n_iter, rows = $n_row, procs = $n_proc, threads = $n_thread)"

	sbatch ${sdir}/run_omp_strong_${N}_${n_iter}_${n_row}_${n_proc}_${n_thread}.slurm
done

#-----------------------------------------------------------------------------#
