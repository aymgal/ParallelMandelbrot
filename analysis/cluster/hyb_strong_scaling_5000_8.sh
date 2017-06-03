#!/bin/bash

#----------------------------- strong scaling --------------------------------#

# initial hybrid strong scaling test

sdir=slurm_runs_strong

N=5000
n_iter=100
n_row=5000  # for best scaling
n_thread=8

outfile=hyb_strong_${N}_${n_iter}_${n_row}_x_${n_thread}.dat
echo "Output file name : " $outfile

mkdir -p outputs

rm -f outputs/$outfile
rm -f out_*.pgm
rm -f out_*.bmp

cd ../../src 			# go there to compile
make hyb
cd ../analysis/cluster	# come back here to run

for n_proc in 2 4 8 12
do
	echo "Going to run several times mandel_hyb, with parameters N = $N, \
max iter = $n_iter, rows = $n_row, procs = $n_proc, threads = $n_thread)"

	sbatch -n ${n_proc} -c ${n_thread} ${sdir}/run_hyb_strong_${N}_${n_iter}_${n_row}_${n_thread}.slurm

done

#-----------------------------------------------------------------------------#
