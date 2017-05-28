#!/bin/bash

rm -rf timings_hyb.dat
rm -rf timings_hyb.png
rm -rf out_*.pgm
rm -rf out_*.bmp

cd ../../src # go there to compile
make hyb
cd ../analysis/shell # return here to analyse

n_nodes=5
n_threads=4
n_iter=100
n_rows=4
for N in 128 256 512 1024 2048 4096 8192
do
	export OMP_NUM_THREADS=$n_threads
    echo "Running mandel_hyb with $n_threads threads on a grid of size $N"
    mpirun -np $n_nodes ../../src/mandel_hyb $N $n_iter $n_rows \
    >> timings_hyb.dat
done

gnuplot << EOF
set terminal png
set output 'timings_hyb.png'
set multiplot layout 1,1
set key off

set logscale x
set logscale y
set xlabel "N [-]"
set ylabel "time [s]"
plot "timings_hyb.dat" using 1:2 with lp
EOF

open timings_hyb.png
