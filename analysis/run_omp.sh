#!/bin/bash

rm -rf timings_omp.dat
rm -rf timings_omp.png
rm -rf out_*.pgm
rm -rf out_*.bmp

cd ../src # go there to compile
make omp
cd ../analysis # return here to analyse

n=4
for N in 128 256 512 1024
do
    export OMP_NUM_THREADS=$n
    echo "Running mandel_omp with $n threads on a grid of size $N"
    ../src/mandel_omp $N 200 >> timings_omp.dat
done

gnuplot <<EOF
set terminal png
set output 'timings_omp.png'
set multiplot layout 1,1
set key off

set logscale y
set xlabel "N [-]"
set ylabel "time [s]"
plot "timings_omp.dat" u 1:2 w lp
EOF

open timings_omp.png
