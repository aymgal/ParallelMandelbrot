#!/bin/bash
#SBATCH --nodes 32
#SBATCH --ntasks 1
#SBATCH --cpus_per_task 1
#SBATCH --account=phpc2017
####SBATCH --reservation=phpc2017

rm -f timetosol_mpi.dat
rm -f out_*.pgm
rm -f out_*.bmp

cd ../src # go there to compile
make mpi
cd ../analysis # return here to analyse

N_size=1000
n_iter=100
n_rows=100
for n_nodes in 1 2 4 8 16 32
do
	srun -N=$n_nodes-$n_nodes mandel_mpi $N_size $n_iter $n_rows \
		>> datafiles/timetosol_mpi.dat
done