#!/bin/bash

#------------------------------ weak scaling ---------------------------------#

sdir=slurm_runs_N

n_iter=200

outfile=srl_N_x_${n_iter}_1_1_1.dat
echo "Output file name : " $outfile

mkdir -p outputs

rm -f outputs/$outfile
rm -f out_*.pgm
rm -f out_*.bmp

cd ../../src 			# go there to compile
make srl
cd ../analysis/cluster	# come back here to run

### problem size doubles when number of processors doubles
### --> N multiplied by sqrt(2)

for N in 1024 1448 2048 2896 4096 5793
do
	echo "Going to run several times mandel_srl, with parameters N = $N, \
max iter = $n_iter, rows = 1, procs = 1, threads = 1)"

	sbatch -n 1 ${sdir}/run_srl_N_${N}_100.slurm

done

#-----------------------------------------------------------------------------#

