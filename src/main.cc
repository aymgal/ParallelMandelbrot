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
#ifdef PARALLEL_MPI
  std::cerr << prog_name << " <grid_size> <n_iter_max> <n_rows>" << std::endl;
#else
  std::cerr << prog_name << " <grid_size> <n_iter_max>" << std::endl;
#endif
  exit(0);
}

static void process_args(int argc, char *argv[], 
                         int& var1, int& var2, int& var3) {
  if (argc < 3)
    usage(argv[0]);

  std::stringstream arg1(argv[1]);
  std::stringstream arg2(argv[2]);
  arg1 >> var1;
  arg2 >> var2;

  if ( (arg1.fail()) || (arg2.fail()) )
    usage(argv[0]);

#ifdef PARALLEL_MPI
  if (argc != 4)
    usage(argv[0]);

  std::stringstream arg3(argv[3]);
  arg3 >> var3;

  if (arg3.fail())
    usage(argv[0]);

#else
  var3 = -1; // whatever we want, won't be used in non-MPI code
#endif 

}


int main(int argc, char *argv[]) {

#ifdef PARALLEL_MPI
  MPI_Init(&argc, &argv);
  int prank, psize;
  MPI_Comm_rank(MPI_COMM_WORLD, &prank);
  MPI_Comm_size(MPI_COMM_WORLD, &psize);
#endif

int N, n_max, n_rows;
process_args(argc, argv, N, n_max, n_rows);

#ifdef PARALLEL_OPENMP
  std::cout << "Number of omp threads : " << omp_get_max_threads()
            << std::endl;
#elif PARALLEL_MPI
  std::cout << "Rank of processor : " << prank << std::endl;
#endif

/* initialization of a mandelbrodt instance */
#ifdef PARALLEL_MPI
  Mandelbrot mandel(N, N, XMIN, XMAX, YMIN, YMAX, n_max, n_rows,
                    MPI_COMM_WORLD);
#else
  Mandelbrot mandel(N, N, XMIN, XMAX, YMIN, YMAX, n_max);
#endif

/* compute Mandelbrot set */
#ifdef OUTPUT_IMAGE
  mandel.run(true);
#else
  mandel.run(false);
#endif
/* ---------------------- */

#ifdef PARALLEL_MPI
  MPI_Finalize();
#endif

  return 0;
}


