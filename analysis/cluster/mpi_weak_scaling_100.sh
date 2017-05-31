#!/bin/bash

#------------------------------ weak scaling ---------------------------------#

sdir=slurm_runs_weak

n_iter=100
n_row=100
n_threads=1

outfile=mpi_weak_x_${n_iter}_${n_row}_x_${n_threads}.dat
echo "Output file name : " $outfile

mkdir -p outputs

rm -f outputs/$outfile
rm -f out_*.pgm
rm -f out_*.bmp

cd ../../src 			# go there to compile
make mpi
cd ../analysis/cluster	# come back here to run

n_proc_min=2
n_proc_max=32
n_proc=$n_proc_min

### problem size doubles when number of processors doubles
### --> N multiplied by sqrt(2)

for N in 1024 1448 2048 2896 4096 5793
do
	echo "Going to run several times mandel_mpi, with parameters N = $N, \
max iter = $n_iter, rows = $n_row, procs = $n_proc, threads = $n_threads)"
    
	sbatch -n $n_proc ${sdir}/run_mpi_weak_${N}_${n_iter}_${n_row}_${n_threads}.slurm

    ((n_proc = n_proc * 2))
done

#-----------------------------------------------------------------------------#

