#include <cmath>
#include <iostream>
#ifdef USE_STDCOMPLEX
#include <complex>
#endif
#ifdef PARALLEL_MPI
#include <vector>
#include <stdexcept>
#endif

#include "mandelbrot.hh"
#include "grid.hh"


#define N_ROWS 30
#define WORK_TAG 2
#define READY_TAG 1
#define END_TAG 0


#if PARALLEL_MPI
Mandelbrot::Mandelbrot(int nx, int ny, 
                       dfloat x_min, dfloat x_max, 
                       dfloat y_min, dfloat y_max,
                       int n_iter, MPI_Comm comm)
: m_global_nx(nx), m_global_ny(ny), 
  m_global_xmin(x_min), m_global_xmax(x_max), 
  m_global_ymin(y_min), m_global_ymax(y_max),
  m_max_iter(n_iter), m_mandel_set(nx, ny),
  m_pdumper(new DumperBinary(m_mandel_set.storage(), comm)),
  m_communicator(comm)
#else
Mandelbrot::Mandelbrot(int nx, int ny, 
                       dfloat x_min, dfloat x_max, 
                       dfloat y_min, dfloat y_max,
                       int n_iter)
: m_global_nx(nx), m_global_ny(ny), 
  m_global_xmin(x_min), m_global_xmax(x_max), 
  m_global_ymin(y_min), m_global_ymax(y_max),
  m_max_iter(n_iter), m_mandel_set(nx, ny),
  m_pdumper(new DumperASCII(m_mandel_set.storage()))
#endif
{
  m_mod_z2_th = 4.0;
  m_value_inside  = 0.0;   // so the set appears black
  m_value_outside = 100.0; 

  m_dx = (m_global_xmax - m_global_xmin) / (dfloat)(m_global_nx - 1);
  m_dy = (m_global_ymax - m_global_ymin) / (dfloat)(m_global_ny - 1);

#ifdef USE_DISTANCE_ESTIM
  m_dist2_th = 1.0e-6;
#endif

#ifdef PARALLEL_MPI
  // get the number of proc and the rank in the proc
  MPI_Comm_rank(m_communicator, &m_prank);
  MPI_Comm_size(m_communicator, &m_psize);

#if defined(MPI_SIMPLE)
  init_allworkers_simple();
#elif defined(MPI_MW_SIMPLE)
  init_master_workers_simple();
#elif defined(MPI_MW_BALANCE)
  {};
#else
#error "MACRO 'MPI_' UNDEFINED"
#endif

#else
  // if non-MPI code, LOCAL and GLOBAL variables are all the same
  m_local_nx = m_global_nx;
  m_local_ny = m_global_ny;
  // and no offsets are needed
  m_local_offset_x = m_local_offset_y = 0;
#endif

}

#if PARALLEL_MPI
void Mandelbrot::init_allworkers_simple() {
  m_MW_communicator = m_communicator;

  // divide the grid on m_psize rows of height m_local_nx
  m_local_nx = m_global_nx / m_psize + (m_prank < m_global_nx % m_psize ? 1 : 0);
  m_local_ny = m_global_ny;
  // offsets
  m_local_offset_x = (m_global_nx / m_psize) * m_prank + 
                     (m_prank < m_global_nx % m_psize ? m_prank : m_global_nx % m_psize);
  m_local_offset_y = 0;

#ifdef VERBOSE
  std::cerr << m_prank << " " 
            << m_global_nx << " " << m_global_ny << " " 
            << m_local_nx << " " << m_local_ny << " " 
            << m_local_offset_x << " " << m_local_offset_y << std::endl;
#endif

  // resizing the grid
  m_mandel_set.resize(m_local_nx, m_local_ny);
}

void Mandelbrot::init_master_workers_simple() {
  // set colors of master and workers communicators
  int color = (prank == 0 ? 0 : 1);
  // create new communicator for master alone or all workers
  MPI_Comm_split(m_communicator, color, m_prank, &m_MW_communicator);

  // define a buffer for local sizes
  std::vector<int> buf_locals(4);

  /*** MASTER ***/
  if (m_prank == 0) {
    int n_workers;
    int w_prank; // worker prank
    int w_nx, w_ny, w_offset_x, w_offset_y;  // worker local problem sizes

    // number of workers
    n_workers = m_psize - 1;
    std::cerr << "Number of workers : " << n_workers << std::endl;

    /* compute local sizes for each worker + put them in lists to send */
    for (int i = 0; i < n_workers; i++) {
      w_prank = i+1;

      // local widths and heights of worker
      w_nx = m_global_nx / n_workers + (i < m_global_nx % n_workers ? 1 : 0);
      w_ny = m_global_ny;
      // local offsets of worker
      w_offset_x = (m_global_nx / n_workers) * i + 
                   (i < m_global_nx % n_workers ? i : m_global_nx % n_workers);
      w_offset_y = 0;

      // put these infos into a single buffer
      buf_locals[0] = w_nx;
      buf_locals[1] = w_ny;
      buf_locals[2] = w_offset_x;
      buf_locals[3] = w_offset_y;

      // and send them to worker
      MPI_Send(&buf_locals[0], 4, MPI_INT, w_prank, WORK_TAG, m_communicator);

      std::cerr << "Sent data to rank " << w_prank << std::endl;

    }
  } /* end master */

  /*** WORKER ***/
  if (m_prank != 0) {
    MPI_Status status;

    // receive buffer from master containing local problem sizes
    MPI_Recv(&buf_locals[0], 4, MPI_INT, 0, WORK_TAG, m_communicator, &status);

    // unpack variables and update member variables
    m_local_nx       = buf_locals[0];
    m_local_ny       = buf_locals[1];
    m_local_offset_x = buf_locals[2];
    m_local_offset_y = buf_locals[3];

    std::cerr << "Rank " << m_prank << " received data from rank 0" 
              << std::endl;
#ifdef VERBOSE
    std::cerr << m_prank << " " 
              << m_global_nx << " " << m_global_ny << " "
              << m_local_nx << " " << m_local_ny << " " 
              << m_local_offset_x << " " << m_local_offset_y << std::endl;
#endif

    // resizing the grid
    m_mandel_set.resize(m_local_nx, m_local_ny);
    
    // set the dumper communicator so that only workers write the image
    m_pdumper->set_communicator(m_MW_communicator);
  } /* end worker */

}


void Mandelbrot::do_all_balance() {
  int worker_ready;
  MPI_Status status;
  // MPI_Request request;

  // set colors of master and workers communicators
  int color = (prank == 0 ? 0 : 1);
  // create new communicator for master alone or all workers
  MPI_Comm_split(m_communicator, color, m_prank, &m_MW_communicator);

  // define a buffer for local sizes
  std::vector<int> buf_locals(4);

  /*** MASTER ***/
  if (m_prank == 0) {
    int i;
    int n_workers;
    int w_prank; // worker prank
    int w_nx, w_ny, w_offset_x, w_offset_y;  // worker local problem sizes
    int row_count;

    // number of workers
    n_workers = m_psize - 1;
    if (n_workers > N_ROWS)
      throw std::invalid_argument("Too much workers processors !");
    std::cerr << "Number of workers : " << n_workers << std::endl;

    /* compute local sizes for each worker + put them in lists to send */
    for (i = 0; i < n_workers; i++) {
      w_prank = i+1;

      // local widths and heights of worker
      w_nx = m_global_nx / N_ROWS + (i < m_global_nx % N_ROWS ? 1 : 0);
      w_ny = m_global_ny;
      // local offsets of worker
      w_offset_x = (m_global_nx / N_ROWS) * i + 
                   (i < m_global_nx % N_ROWS ? i : m_global_nx % N_ROWS);
      w_offset_y = 0;

      // put these infos into a single buffer
      buf_locals[0] = w_nx;
      buf_locals[1] = w_ny;
      buf_locals[2] = w_offset_x;
      buf_locals[3] = w_offset_y;

      // and send them to worker
      MPI_Send(&buf_locals[0], 4, MPI_INT, w_prank, WORK_TAG, m_communicator);
      std::cerr << "(init) Rank 0 send data to rank " << w_prank << std::endl;
    }

    /* 'infinite' loop to feed workers with new rows */
    n_busy_workers = n_workers;
    row_count = n_workers; // number of already computed/written rows
    for (;;) {
      worker_ready = 0;
      MPI_Recv(&worker_ready, 1, MPI_INT, 
               MPI_ANY_SOURCE, READY_TAG, m_communicator, &status);

      w_prank = status.MPI_SOURCE;
      std::cerr << "Rank 0 recv 'ready' from " << w_prank << std::endl;

      if ((worker_ready == 1) && (row_count < N_ROWS)) {
        n_busy_workers--;
        row_count++;

        row_idx = row_count - 1; // index of current row

        // local widths and heights of worker
        w_nx = m_global_nx / N_ROWS + (row_idx < m_global_nx % N_ROWS ? 1 : 0);
        w_ny = m_global_ny;
        // local offsets of worker
        w_offset_x = (m_global_nx / N_ROWS) * row_idx + 
                     (row_idx < m_global_nx % N_ROWS ? row_idx : m_global_nx % N_ROWS);
        w_offset_y = 0;

        // put these infos into a single buffer
        buf_locals[0] = w_nx;
        buf_locals[1] = w_ny;
        buf_locals[2] = w_offset_x;
        buf_locals[3] = w_offset_y;

        // and send them to worker
        MPI_Send(&buf_locals[0], 4, MPI_INT, w_prank, WORK_TAG, m_communicator);
        std::cerr << "Rank 0 send data to rank " << w_prank << std::endl;
      }

      if (row_count == N_ROWS) {
        std::cerr << "ALL ROWS COMPLETE" << std::endl;
        break;
      }
    }

    // when all rows are complete, send a message to all workers to quit
    for (i = 0; i < n_workers; i++) {
      w_prank = i+1;
      buf_locals = {0, 0, 0, 0}; // whatever we want
      MPI_Send(&buf_locals[0], 4, MPI_INT, w_prank, END_TAG, m_communicator);
    }
  
  } /* end master */


  /*** WORKER ***/
  if (m_prank != 0) {

    // set the dumper communicator so that only workers write the image
    m_pdumper->set_communicator(m_MW_communicator);

    for (;;) {
      // receive buffer from master containing local problem sizes
      MPI_Recv(&buf_locals[0], 4, MPI_INT, 
               0, MPI_ANY_TAG, m_communicator, &status);

      if (status.MPI_TAG == END_TAG) {
        std::cerr << "END TAG RECV" << std::endl;
        return;
      }

      worker_ready = 0;

      // unpack variables and update member variables
      m_local_nx       = buf_locals[0];
      m_local_ny       = buf_locals[1];
      m_local_offset_x = buf_locals[2];
      m_local_offset_y = buf_locals[3];

      std::cerr << "Rank " << m_prank << " recv data from rank 0" 
                << std::endl;
  #ifdef VERBOSE
      std::cerr << m_prank << " " 
                << m_global_nx << " " << m_global_ny << " "
                << m_local_nx << " " << m_local_ny << " " 
                << m_local_offset_x << " " << m_local_offset_y << std::endl;
  #endif

      // resizing the grid
      m_mandel_set.resize(m_local_nx, m_local_ny);

      /* compute the set and write corresponding part in image file */
      compute_set();
      write_image(m_local_offset_x, m_global_nx);

      /* finally send message to tell master that worker is ready for new work */
      worker_ready = 1;
      MPI_Send(&worker_ready, 1, MPI_INT, 0, READY_TAG, m_communicator);
      std::cerr << "Rank " << m_prank << " send 'ready' to rank 0" 
                << std::endl;
    }

  } /* end worker */
}

#endif /* PARALLEL_MPI */


void Mandelbrot::compute_set() {
  int ix, iy;

#ifdef MPI_MW_SIMPLE
  if (m_prank == 0) return; // master does not compute
#endif

#ifdef PARALLEL_OPENMP
  /* OpenMP parallelization is set here */
  #pragma omp parallel for private(ix, iy) schedule(dynamic)
#endif

  for (ix = 0; ix < m_local_nx; ix++) {

    for (iy = 0; iy < m_local_ny; iy++) {

      /* compute current pixel */
      compute_pix(ix, iy);

    }
  }

}


void Mandelbrot::compute_pix(int ix, int iy) {
  dfloat cx, cy;
  dfloat z0x, z0y;

  /* compute complex coordinates of current point */
  cx = m_global_xmin + (dfloat)(ix + m_local_offset_x) * m_dx;
  cy = m_global_ymin + (dfloat)(iy + m_local_offset_y) * m_dy;

  /* set initial condition of iteration */
  z0x = z0y = 0.0;

  Grid & mset = m_mandel_set.storage();
  /* call to the main iteration loop on current pixel */
  mset(ix, iy) = iterate_on_pix(cx, cy, z0x, z0y);

}

dfloat Mandelbrot::iterate_on_pix(dfloat cx, dfloat cy, 
                                  dfloat z0x, dfloat z0y) {
  int iter;
  dfloat zx, zy;
  dfloat mod_z2;
  bool diverge;
#ifdef USE_STDCOMPLEX
  std::complex<dfloat> c;
  std::complex<dfloat> z, z2;
  std::complex<dfloat> z_new;
#else
  dfloat zx2, zy2;
  dfloat zx_new, zy_new;
#endif
#ifdef USE_DISTANCE_ESTIM
  dfloat dzx, dzy, dzx_new, dzy_new;
  dfloat mod_dz2;
  dfloat dist2;
#endif
  dfloat pix_value;

/* initialization before iteration loop */
#ifdef USE_STDCOMPLEX
  c = {cx, cy};
  z = {z0x, z0y};
#else
  zx = z0x;
  zy = z0y;
#endif
  mod_z2 = 0.0;
  diverge = false;
#ifdef USE_DISTANCE_ESTIM
  dzx = dzy = dzx_new = dzy_new = 0.0;
  mod_dz2 = 0.0;
  dist2 = 0.0;
#endif

  /*** MAIN LOOP that computes if a point belong to the set or not ***/
  for (iter = 1; iter <= m_max_iter; iter++) {

    /* compute the recursive equation */
#ifdef USE_STDCOMPLEX
    // trivial if using std::complex
    z_new = z*z + c;

    // less trivial otherwise...
#else
    // save squared coordinates to be reused
    zx2 = zx*zx;
    zy2 = zy*zy;
    // recursive equation, real part
    zx_new = zx2 - zy2 + cx;
    // recursive equation, imaginary part
    // but zy_new = 2 * zx * zy + cy computed in a 'clever' way
#ifdef SQUARE_TRICK // <-- this should be efficient only at large scales
    zy_new = pow(zx * zy, 2) - zx2 - zy2;
#else
    zy_new = zx * zy;
    zy_new += zy_new; // replaces multiplication by 2
#endif
    zy_new += cy;
#endif

#ifdef USE_DISTANCE_ESTIM
    /* compute derivative for distance estimator */
    // in a 'clever' way, in place of
    // dzx_new = 2*(zx*dzx-zy*dzy) + 1 and dzy_new = 2*(zx*dzy+zy*dzx)
#ifdef USE_STDCOMPLEX
    zx = z.real();
    zy = z.imag();
#endif
    dzx_new = zx*dzx - zy*dzy;
    dzx_new += dzx_new;
    dzx_new += 1.0;
    dzy_new = zx*dzy + zy*dzx;
    dzy_new += dzy_new;
#endif

    /* update variables */
#ifdef USE_STDCOMPLEX
    z = z_new;
#else
    zx = zx_new;
    zy = zy_new;
#endif

#ifdef USE_DISTANCE_ESTIM
    dzx = dzx_new;
    dzy = dzy_new;
#endif

    // compute modulus
#ifdef USE_STDCOMPLEX
    mod_z2 = std::norm(z);
#else
    mod_z2 = zx2 + zy2;
#endif

    /* check if point diverges or not */
    if (mod_z2 > m_mod_z2_th) {

      diverge = true;

#ifdef USE_DISTANCE_ESTIM
      mod_dz2 = dzx*dzx + dzy*dzy;
      dist2   = pow(log(mod_z2), 2) * mod_z2 / mod_dz2;
#endif
      break;
    }

  } /*** end main loop ***/


  /* now, assign value to current pixel depending on some conditions */

  // point 'inside' the set
#ifdef USE_DISTANCE_ESTIM
  if ((!diverge) || (dist2 < m_dist2_th)) {
    pix_value = m_value_inside;
  }
#else
  if (!diverge) {
    pix_value = m_value_inside;
  }
#endif

  // point 'outside' the set
  else {
#ifdef COLOR_ITERATIONS
    pix_value = (dfloat) iter;
#elif defined(COLOR_PRANK) && defined(PARALLEL_MPI)
    pix_value = m_prank + 1;
#else
    pix_value = m_value_outside;
#endif
  }

  return pix_value;
}


void Mandelbrot::write_image(int arg1, int arg2) {
#if defined(PARALLEL_MPI) && (defined(MPI_MW_SIMPLE) || defined(MPI_MW_BALANCE))
  if (m_prank == 0) return;
#endif
  m_pdumper->dump(arg1, arg2);
}


