#ifndef MANDEL_HH
#define MANDEL_HH

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
             int n_iter, MPI_Comm comm);
#else
  Mandelbrot(int nx, int ny, 
           dfloat x_min, dfloat x_max, 
           dfloat y_min, dfloat y_max,
           int n_iter);
#endif

  void compute_set();

  void write_image(int arg1, int arg2);

  void do_all_balance();

private:
  // global size of the problem
  int m_global_nx, m_global_ny;
  dfloat m_global_xmin, m_global_xmax;
  dfloat m_global_ymin, m_global_ymax;

  // max number of iterations
  int m_max_iter;

  // grid storage
  Buffer m_mandel_set;

  // dumper to use for outputs
  std::unique_ptr<Dumper> m_pdumper;

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


#ifdef PARALLEL_MPI

  // proc rank
  int m_prank;
  // communicator size
  int m_psize;

  // communicators
  MPI_Comm m_communicator;    // main communicator
  MPI_Comm m_MW_communicator; // master/workers communicator

  /* '_simple' : space divided in equal psize rows */
  void init_allworkers_simple();
  void init_master_workers_simple();

  /* '_balance' : space divided in a lot of rows, 
   * then master ditribute rows to workers dynamically
   */
   // void master_balance();
   // void worker_balance();
  // void init_master_workers_balance();

  /* for workers, to compute then write */
  // void compute_and_write();

#endif /* PARALLEL_MPI */

  void compute_pix(int ix, int iy);
  dfloat iterate_on_pix(dfloat cx, dfloat cy, dfloat z0x, dfloat z0y);

};


#endif /* MANDEL_HH */

