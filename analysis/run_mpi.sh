#!/bin/bash

rm -rf timings_mpi.dat
rm -rf timings_mpi.png
rm -rf out_*.pgm
rm -rf out_*.bmp

cd ../src # go there to compile
# make mpi
cd ../analysis # return here to analyse

n_nodes=5
n_iter=200
n_rows=4
output=1
for N in 128 256 512 1024 2048 4096
do
    echo "Running mandel_mpi on a grid of size $N"
    mpirun -np $n_nodes ../src/mandel_mpi $N $n_iter $n_rows $output \
    	>> timings_mpi.dat
done

gnuplot <<EOF
set terminal png
set output 'timings_mpi.png'
set multiplot layout 1,1
set key off

set logscale x
set logscale y
set xlabel "N [-]"
set ylabel "time [s]"
plot "timings_mpi.dat" u 1:2 w lp
EOF

open timings_mpi.png
