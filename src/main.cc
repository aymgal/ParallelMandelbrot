#include <iostream>
#include <sstream>

#include "gvars.hh"
#include "mandelbrot.hh"

#ifdef PARALLEL_OPENMP
#include <omp.h> // for omp_get_max_threads()
#elif PARALLEL_MPI
#include <mpi.h>
#endif


static void usage(const std::string & prog_name) {
  std::cerr << prog_name << " <grid_size> <n_iter_max>" << std::endl;
  exit(0);
}

int main(int argc, char *argv[]) {

#ifdef PARALLEL_MPI
  MPI_Init(&argc, &argv);
  int prank, psize;
  MPI_Comm_rank(MPI_COMM_WORLD, &prank);
  MPI_Comm_size(MPI_COMM_WORLD, &psize);
#endif

  if (argc != 3) usage(argv[0]);

  std::stringstream arg1(argv[1]);
  std::stringstream arg2(argv[2]);
  int N, n_max; // main two variables of the problem
  arg1 >> N;
  arg2 >> n_max;

  if ((arg1.fail()) || (arg2.fail())) usage(argv[0]);

#ifdef PARALLEL_OPENMP
  std::cout << "Number of omp threads : " << omp_get_max_threads()
            << std::endl;
#elif PARALLEL_MPI
  std::cout << "Rank of processor : " << prank << std::endl;
#endif

/* initialization of a mandelbrodt instance */
#ifdef PARALLEL_MPI
  Mandelbrot mandel(N, N, XMIN, XMAX, YMIN, YMAX, n_max, MPI_COMM_WORLD);
#else
  Mandelbrot mandel(N, N, XMIN, XMAX, YMIN, YMAX, n_max);
#endif

/* start timings */
#ifdef PARALLEL_MPI
  MPI_Barrier(MPI_COMM_WORLD);
  auto mpi_start_c = MPI_Wtime();
#endif
  auto start_c = clk::now();

/* compute Mandelbrot set */
#ifdef MPI_MW_BALANCE
  mandel.do_all_balance();
#else
  mandel.compute_set();
#endif
/* ---------------------- */

/* end timings */
  auto end_c = clk::now();
#ifdef PARALLEL_MPI
  MPI_Barrier(MPI_COMM_WORLD);
  auto mpi_end_c = MPI_Wtime();
#endif

  second time_compute = end_c - start_c;
  std::cout << "Time to compute set : ";
#ifdef PARALLEL_MPI
  std::cout << prank << " " << N << " " << time_compute.count() << std::endl;
#else
  std::cout << N << " " << time_compute.count() << std::endl;
#endif

#ifdef PARALLEL_MPI
  if (prank == 0) 
    std::cout << "MPI time to compute " << mpi_end_c-mpi_start_c << std::endl; 
#endif

#if defined(OUTPUT_IMAGE) && !defined(MPI_MW_BALANCE)

  auto start_w = clk::now();

/* output Mandelbrot set as an image */
  mandel.write_image(0);
/* --------------------------------- */

  auto end_w = clk::now();

  second time_write = end_w - start_w;
  std::cout << "Time to write image : ";
#ifdef PARALLEL_MPI
  std::cout << prank << " " << N << " " << time_write.count() << std::endl;
#else
  std::cout << N << " " << time_write.count() << std::endl;
#endif

#endif

#ifdef PARALLEL_MPI
  MPI_Finalize();
#endif

  return 0;
}
