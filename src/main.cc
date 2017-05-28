#include <iostream>
#include <sstream>

#include "gvars.hh"
#include "mandelbrot.hh"

#ifdef PARALLEL_OPENMP
#include <omp.h> // for omp_get_max_threads()
#endif
#if PARALLEL_MPI
#include <mpi.h>
#endif


static void usage(const std::string & prog_name);

#ifdef PARALLEL_MPI
static void process_args(int argc, char *argv[], 
                         int& int1, int& int2, int& int3, bool& bool1);
#else
static void process_args(int argc, char *argv[], 
                         int& int1, int& int2, bool& bool1);
#endif


int main(int argc, char *argv[]) {

#ifdef PARALLEL_MPI

#ifdef PARALLEL_OPENMP // for hybrid program (= multi-threading MPI)
  int required, provided;
  required = MPI_THREAD_FUNNELED;
  MPI_Init_thread(&argc, &argv, required, &provided);
#else
  MPI_Init(&argc, &argv);
#endif

  int prank, psize;
  MPI_Comm_rank(MPI_COMM_WORLD, &prank);
  MPI_Comm_size(MPI_COMM_WORLD, &psize);
#endif

#if defined(VERBOSE) && defined(PARALLEL_OPENMP)
  std::cout << "Number of omp threads : " << omp_get_max_threads()
            << std::endl;
#elif defined(VERBOSE) && defined(PARALLEL_MPI)
  std::cout << "Rank of processor : " << prank << std::endl;
#elif defined(VERBOSE)
  std::cout << "Serial code" << prank << std::endl;
#endif

/* initialization of a mandelbrodt instance */
#ifdef PARALLEL_MPI
  int N, n_max, n_rows;
  bool output_img;
  process_args(argc, argv, N, n_max, n_rows, output_img);

  MandelbrotSet mandel(N, N, XMIN, XMAX, YMIN, YMAX, n_max, output_img, 
                       n_rows, MPI_COMM_WORLD);
#else
  int N, n_max;
  bool output_img;
  process_args(argc, argv, N, n_max, output_img);

  MandelbrotSet mandel(N, N, XMIN, XMAX, YMIN, YMAX, n_max, output_img);
#endif

  /* compute Mandelbrot set and write an image if set so */
  mandel.run();

#ifdef PARALLEL_MPI
  MPI_Finalize();
#endif

  return 0;
}


static void usage(const std::string & prog_name) {
#ifdef PARALLEL_MPI
  std::cerr << prog_name << " <grid_size> <n_iter_max> <n_rows> <output_img>"
            << std::endl;
#else
  std::cerr << prog_name << " <grid_size> <n_iter_max> <output_img>" 
           << std::endl;
#endif
  exit(0);
}

#ifdef PARALLEL_MPI
static void process_args(int argc, char *argv[], 
                          int& int1, int& int2, int& int3, bool& bool1) {
  if (argc < 5)
    usage(argv[0]);

  std::stringstream arg1(argv[1]);
  std::stringstream arg2(argv[2]);
  std::stringstream arg3(argv[3]);
  std::stringstream arg4(argv[4]);
  arg1 >> int1;
  arg2 >> int2;
  arg3 >> int3;
  arg4 >> bool1;

  if ( (arg1.fail()) || (arg2.fail()) || (arg3.fail()) || (arg4.fail()) )
    usage(argv[0]);
}

#else

static void process_args(int argc, char *argv[], 
                         int& int1, int& int2, bool& bool1) {
  if (argc < 4)
    usage(argv[0]);

  std::stringstream arg1(argv[1]);
  std::stringstream arg2(argv[2]);
  std::stringstream arg3(argv[3]);
  arg1 >> int1;
  arg2 >> int2;
  arg3 >> bool1;

  if ( (arg1.fail()) || (arg2.fail()) || (arg3.fail()) )
    usage(argv[0]);
}
#endif
