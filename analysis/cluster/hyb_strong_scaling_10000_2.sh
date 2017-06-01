#!/bin/bash

#----------------------------- strong scaling --------------------------------#

# initial hybrid strong scaling test

sdir=slurm_runs_strong

N=10000
n_iter=100
n_row=1  # because only 1 worker (because n_proc=2)
n_proc=2

outfile=hyb_strong_${N}_${n_iter}_${n_row}_${n_proc}_x.dat
echo "Output file name : " $outfile

mkdir -p outputs

rm -f outputs/$outfile
rm -f out_*.pgm
rm -f out_*.bmp

cd ../../src 			# go there to compile
make hyb
cd ../analysis/cluster	# come back here to run

for n_thread in 1 2 4 8 16
do
	echo "Going to run several times mandel_hyb, with parameters N = $N, \
max iter = $n_iter, rows = $n_row, procs = $n_proc, threads = $n_thread)"
	
	if [[ ${n_thread} -le 4 ]]
	then
		sbatch -n ${n_proc} ${sdir}/run_hyb_strong_${N}_${n_iter}_${n_row}_${n_proc}_${n_thread}.slurm
	else
		sbatch -n ${n_proc} -c ${n_thread} \
			${sdir}/run_hyb_strong_${N}_${n_iter}_${n_row}_${n_proc}_x.slurm
	fi

done

#-----------------------------------------------------------------------------#
