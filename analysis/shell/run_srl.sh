#!/bin/bash

rm -rf timings_srl.dat
rm -rf timings_srl.png
rm -rf out_*.pgm
rm -rf out_*.bmp

cd ../../src # go there to compile
make srl
cd ../analysis/shell # return here to analyse

n_iter=100
for N in 128 256 512 1024 2048 4096 8192
do
    echo "Running mandel_srl on a grid of size $N"
    ../../src/mandel_srl $N $n_iter >> timings_srl.dat
done

gnuplot << EOF
set terminal png
set output 'timings_srl.png'
set multiplot layout 1,1
set key off

set logscale x
set logscale y
set xlabel "N [-]"
set ylabel "time [s]"
plot "timings_srl.dat" using 1:2 with lp
EOF

open timings_srl.png
