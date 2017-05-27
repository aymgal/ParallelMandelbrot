#ifndef MANDEL_HH
#define MANDEL_HH

#include <sstream>
#include <vector>

#include "gvars.hh"
#include "buffer.hh"
#include "dumpers.hh"

#ifdef PARALLEL_MPI
#include <mpi.h>
#endif

class MandelbrotSet {
public:
#ifdef PARALLEL_MPI
  MandelbrotSet(int nx, int ny, 
                dfloat x_min, dfloat x_max, 
                dfloat y_min, dfloat y_max,
                int n_iter, bool output_img, 
                int n_rows, MPI_Comm comm);
#else
  MandelbrotSet(int nx, int ny, 
                dfloat x_min, dfloat x_max, 
                dfloat y_min, dfloat y_max,
                int n_iter, bool output_img);
#endif

  void run();

private:
  // global size of the problem
  int m_global_nx, m_global_ny;
  dfloat m_global_xmin, m_global_xmax;
  dfloat m_global_ymin, m_global_ymax;

  // max number of iterations
  int m_max_iter;

  // toggle on/off output image
  bool m_output_img;

  // grid storage
  Buffer m_mandel_set;

  // dumper to use for outputs
#ifdef PARALLEL_MPI
  std::unique_ptr<DumperBinary> m_pdumper;
#else
  std::unique_ptr<DumperASCII> m_pdumper;
#endif

  // threshold squared modulus
  dfloat m_mod_z2_th;

  // values of pixels belonging (inside) or not (outside) to the set
  dfloat m_value_inside, m_value_outside;

  // grid step in each directions
  dfloat m_dx, m_dy;

#ifdef USE_DISTANCE_ESTIM
  // threshold distance
  dfloat m_dist2_th;
#endif

  // local (on current proc) problem size
  int m_local_nx, m_local_ny;
  
  // local offset
  int m_local_offset_x, m_local_offset_y;

  // compute the Mandelbrot set
  void compute_set();

  // assign value to pixel in the grid
  void compute_pix(int ix, int iy);

  // compute the complex recursive equation, returns value of the pixel
  dfloat solve_recursive(dfloat cx, dfloat cy, dfloat z0x, dfloat z0y);

  // compute local sizes/offsets defining a row for the grid decomposition
  std::vector<int> get_row_def(int row_idx, int nx, int ny, int n_rows);

  void init_dumper_colors();

#ifdef PARALLEL_MPI
  // proc rank
  int m_prank;

  // communicator size
  int m_psize;

  // number of rows to divide complex grid
  int m_n_rows;

  // communicators
  MPI_Comm m_communicator;    // main communicator
  MPI_Comm m_MW_communicator; // (master) and (all workers) communicator

  // initialize local sizes so every proc compute a given row
  void init_mpi_simple();

  // for master/workers load balancing
  void init_writers(int prank_nonwriter);
  void mpi_master();
  void mpi_worker();

#endif /* PARALLEL_MPI */

};


#endif /* MANDEL_HH */

