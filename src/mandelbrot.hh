#ifndef MANDEL_HH
#define MANDEL_HH

#include <vector>

#include "gvars.hh"
#include "buffer.hh"
#include "dumpers.hh"

#ifdef PARALLEL_MPI
#include <mpi.h>
#endif

class Mandelbrot {
public:
#ifdef PARALLEL_MPI
  Mandelbrot(int nx, int ny, 
             dfloat x_min, dfloat x_max, 
             dfloat y_min, dfloat y_max,
             int n_iter, int n_rows, MPI_Comm comm);
#else
  Mandelbrot(int nx, int ny, 
             dfloat x_min, dfloat x_max, 
             dfloat y_min, dfloat y_max,
             int n_iter);
#endif

  void run(bool output_img);

private:
  // global size of the problem
  int m_global_nx, m_global_ny;
  dfloat m_global_xmin, m_global_xmax;
  dfloat m_global_ymin, m_global_ymax;

  // max number of iterations
  int m_max_iter;

  // grid storage
  Buffer m_mandel_set;

#ifdef PARALLEL_MPI
  // dumper to use for outputs
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

  std::vector<int> get_row_def(int row_idx, int nx, int ny, int n_rows);

#ifdef PARALLEL_MPI
  // proc rank
  int m_prank;
  // communicator size
  int m_psize;

  // number of rows to divide complex grid
  int m_n_rows;

  // communicators
  MPI_Comm m_communicator;    // main communicator
  MPI_Comm m_MW_communicator; // master/workers communicator

  /* '_simple' : complex plane divided in equal psize rows */
  void init_mpi_simple();      // to be called in constructor

  /* '_balance' : complex plane divided in a lot of rows, 
   * then master ditribute rows to workers dynamically
   */
   void mpi_master();
   void mpi_worker(bool output_img);

#endif /* PARALLEL_MPI */

};


#endif /* MANDEL_HH */

