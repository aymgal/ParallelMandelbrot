#!/bin/bash

#------------------------------ weak scaling ---------------------------------#

n_iter=100
n_row=100
n_threads=0	# because unused in the case of non-OpenMP code

outfile=mpi_weak_x_${n_iter}_${n_row}_x_${n_threads}.dat
echo "Output file name : " $outfile

rm -f outputs/$outfile
rm -f out_*.pgm
rm -f out_*.bmp

cd ../../src 			# go there to compile
# make mpi
cd ../analysis/cluster	# come back here to run

N_min=512
N_max=8192
N=$N_min
n_proc_min=2
n_proc_max=32
n_proc=$n_proc_min

while [[ $N -le $N_max ]] && [[ $n_proc -le $n_proc_max ]]
do
	echo "Going to run several times mandel_mpi, with parameters N = $N, \
max iter = $n_iter, rows = $n_row, procs = $n_proc, threads = $n_threads)"
    
	sbatch -n $n_proc run_mpi_soft_N${N}.slurm

    ((N = N * 2))			# grows in the same time system size...
    ((n_proc = n_proc * 2))	# ... and number of processes
done

#-----------------------------------------------------------------------------#

