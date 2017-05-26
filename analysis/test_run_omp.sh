#! /bin/bash

for n in {2..4}
do
    export OMP_NUM_THREADS=$n
    echo "Running mandel_omp with $n threads"
    ../src/mandel_omp 1000 200
done