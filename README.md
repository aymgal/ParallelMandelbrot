# Parallelized Mandelbrot set generator

**Parallel implementations of Mandelbrot set program, written in a C++ object-oriented approach.**

Different versions :
* serial
* OpenMP
* MPI
* MPI with master/workers load-balancing
* hybrid OpenMP + MPI (not ready yet)

It is written in unified source code. The makefile controls which version is compiled :
* `make srl` for serial version
* `make omp` for OpenMP version
* `make mpi` with flag `MPI_SIMPLE` for simple MPI version
* `make mpi` with flag `MPI_MASTER_WORKERS` for load-balancing MPI
* `make` for all.
