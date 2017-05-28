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


// define limits in the complex plane
#define XMIN -2.0
#define XMAX  1.0
#define YMIN -1.5
#define YMAX  1.5


static void usage(const std::string & prog_name);

#if defined(PARALLEL_MPI) && defined(MPI_MASTER_WORKERS)
static void process_args(int argc, char *argv[], 
                         int& int1, int& int2, int& int3);
#else
static void process_args(int argc, char *argv[], 
                         int& int1, int& int2);
#endif

/* ------------------------------------------------------------------------- */

int main(int argc, char *argv[]) {

#ifdef PARALLEL_MPI
#ifdef PARALLEL_OPENMP // for hybrid program (= multi-threading MPI)
  int required, provided;
  required = MPI_THREAD_FUNNELED;
  // required = MPI_THREAD_SERIALIZED;
  // required = MPI_THREAD_MULTIPLE;
  MPI_Init_thread(&argc, &argv, required, &provided);
#ifdef VERBOSE
    std::cerr << "Levels of thread support " 
              << required << " " << provided << std::endl;
#endif
#else
  MPI_Init(&argc, &argv);
#endif

  int prank, psize;
  MPI_Comm_rank(MPI_COMM_WORLD, &prank);
  MPI_Comm_size(MPI_COMM_WORLD, &psize);
#endif

#ifdef VERBOSE
#if defined(PARALLEL_OPENMP) && defined(PARALLEL_MPI)
  std::cout << "Number of omp threads for rank " << prank << " : "
            << omp_get_max_threads() << std::endl;
#elif defined(PARALLEL_OPENMP)
  std::cout << "Number of omp threads : " << omp_get_max_threads()
            << std::endl;
#elif && defined(PARALLEL_MPI)
  std::cout << "Rank of processor : " << prank << std::endl;
#else
  std::cout << "Serial code" << prank << std::endl;
#endif
#endif

/* initialization of a mandelbrodt instance */
  int N, n_max;
#if defined(PARALLEL_MPI) && defined(MPI_MASTER_WORKERS)
  int n_rows;
  process_args(argc, argv, N, n_max, n_rows);

  MandelbrotSet mandel(N, N, XMIN, XMAX, YMIN, YMAX, n_max, 
                       n_rows, MPI_COMM_WORLD);
  
#elif defined(PARALLEL_MPI) && defined(MPI_SIMPLE)
  int n_rows = 0; // useless in this case (non-load balanced MPI program)
  process_args(argc, argv, N, n_max);

  MandelbrotSet mandel(N, N, XMIN, XMAX, YMIN, YMAX, n_max, 
                       n_rows, MPI_COMM_WORLD);

#else
  process_args(argc, argv, N, n_max);

  MandelbrotSet mandel(N, N, XMIN, XMAX, YMIN, YMAX, n_max);
#endif

  /* compute Mandelbrot set and write an image if set so */
  mandel.run();

#ifdef PARALLEL_MPI
  MPI_Finalize();
#endif

  return 0;
}

/* ------------------------------------------------------------------------- */

static void usage(const std::string & prog_name) {
#ifdef PARALLEL_MPI
  std::cerr << prog_name << " <grid_size> <n_iter_max> <n_rows>"
            << std::endl;
#else
  std::cerr << prog_name << " <grid_size> <n_iter_max>" 
           << std::endl;
#endif
  exit(0);
}

#if defined(PARALLEL_MPI) && defined(MPI_MASTER_WORKERS)
static void process_args(int argc, char *argv[], 
                         int& int1, int& int2, int& int3) {
  if (argc < 4)
    usage(argv[0]);

  std::stringstream arg1(argv[1]);
  std::stringstream arg2(argv[2]);
  std::stringstream arg3(argv[3]);
  arg1 >> int1;
  arg2 >> int2;
  arg3 >> int3;

  if ( (arg1.fail()) || (arg2.fail()) || (arg3.fail()) )
    usage(argv[0]);
}

#else

static void process_args(int argc, char *argv[], 
                         int& int1, int& int2) {
  if (argc < 3)
    usage(argv[0]);

  std::stringstream arg1(argv[1]);
  std::stringstream arg2(argv[2]);
  arg1 >> int1;
  arg2 >> int2;

  if ( (arg1.fail()) || (arg2.fail()) )
    usage(argv[0]);
}
#endif
