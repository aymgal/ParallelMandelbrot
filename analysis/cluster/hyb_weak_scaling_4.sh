#!/bin/bash

#------------------------------ weak scaling ---------------------------------#

# vary number of rows of the division, to be equal to N

sdir=slurm_runs_weak

n_iter=100
n_thread=4 # per rank

outfile=hyb_weak_x_${n_iter}_x_x_${n_thread}.dat
echo "Output file name : " $outfile

mkdir -p outputs

rm -f outputs/$outfile
rm -f out_*.pgm
rm -f out_*.bmp

cd ../../src 			# go there to compile
make hyb
cd ../analysis/cluster	# come back here to run

n_proc_min=2
n_proc_max=32
n_proc=$n_proc_min

### problem size doubles when number of processors doubles
### --> N multiplied by sqrt(2)

for N in 1024 1448 2048 2896 4096 5793
do
	n_row=$N # division in rows equal to N (total of N rows of 1 pixel thick)

	echo "Going to run several times mandel_hyb, with parameters N = $N, \
max iter = $n_iter, rows = $n_row, procs = $n_proc, threads = $n_thread)"
    
	sbatch -n $n_proc -c $n_thread ${sdir}/run_hyb_weak_${N}_${n_iter}_${n_row}_${n_thread}.slurm

	((n_proc = n_proc * 2))
done

#-----------------------------------------------------------------------------#

