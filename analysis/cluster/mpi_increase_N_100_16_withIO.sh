#!/bin/bash

#------------------------------ weak scaling ---------------------------------#

sdir=slurm_runs_N

n_iter=100
n_proc=16
n_thread=1

outfile=mpi_N_x_${n_iter}_x_${n_proc}_1_IO.dat
echo "Output file name : " $outfile

mkdir -p outputs

rm -f outputs/$outfile
rm -f out_*.pgm
rm -f out_*.bmp

cd ../../src 			# go there to compile
make mpi
cd ../analysis/cluster	# come back here to run

### problem size doubles when number of processors doubles
### --> N multiplied by sqrt(2)

for N in 1024 1448 2048 2896 4096
do
	n_row=$N

	echo "Going to run several times mandel_mpi with I/O, \
with parameters N = $N, max iter = $n_iter, rows = $n_row, procs = $n_proc, threads = $n_thread)"

	sbatch -n ${n_proc} ${sdir}/run_mpi_N_${N}_${n_iter}_${n_row}_${n_proc}_IO.slurm

done

#-----------------------------------------------------------------------------#

