# Parallelized Mandelbrot set generator

**Parallel implementations of Mandelbrot set program, written in a C++ object-oriented approach.**

Different versions :
* serial
* OpenMP
* MPI
* MPI with master/workers load-balancing
* hybrid OpenMP + MPI

Output example, for a window in the complex plane defined by [-0.2, 0.4] x [0.5, 1.1] :
<p align="center">
<img src="images/out_4000_100.bmp" alt="Output example" width="400" align="center"/>
</p>

It is written in a unique source code. The Makefile controls which version is compiled :
* `make srl` for serial version
* `make omp` for OpenMP version
* `make mpi` 
	* with `MPI_SIMPLE` flag for simple MPI version
	* with `MPI_MASTER_WORKERS` flag for load-balancing MPI version
* `make hyb` with appropriate MPI flag for hybrid OpenMP + MPI version
* `make` for all.
