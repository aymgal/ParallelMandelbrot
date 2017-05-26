#include <cmath>
#include <iostream>
#include <stdexcept>
#ifdef USE_STDCOMPLEX
#include <complex>
#endif

#include "mandelbrot.hh"
#include "grid.hh"
#include "utils.hh"


#ifdef MPI_MASTER_WORKERS
#define FIRST_WORK_TAG 4
#define WORK_TAG       3
#define FEEBACK_TAG    2
#define END_TAG        1
#define LAST_END_TAG   0
#endif


#if PARALLEL_MPI
Mandelbrot::Mandelbrot(int nx, int ny, 
                       dfloat x_min, dfloat x_max, 
                       dfloat y_min, dfloat y_max,
                       int n_iter, int n_rows, MPI_Comm comm)
    : m_global_nx(nx), m_global_ny(ny), 
      m_global_xmin(x_min), m_global_xmax(x_max), 
      m_global_ymin(y_min), m_global_ymax(y_max),
      m_max_iter(n_iter), m_mandel_set(nx, ny),
      m_pdumper(new DumperBinary(m_mandel_set.storage(), comm)),
      m_n_rows(n_rows), m_communicator(comm)
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
  init_mpi_simple();
#elif defined(MPI_MASTER_WORKERS)
  {};
#else
#error "MACRO 'MPI_' UNDEFINED"
#endif

#else
  // if non-MPI code, LOCAL and GLOBAL variables are the same
  m_local_nx = m_global_nx;
  m_local_ny = m_global_ny;
  // and no offsets are needed
  m_local_offset_x = m_local_offset_y = 0;
#endif

}


void Mandelbrot::run(bool output_img) {

#if defined(PARALLEL_MPI) && defined(MPI_MASTER_WORKERS)
  // create new communicators for (master alone) and (all workers)
  int color = (m_prank == 0 ? 0 : 1);
  MPI_Comm_split(m_communicator, color, m_prank, &m_MW_communicator);

  /*** MPI master ***/
  if (m_prank == 0) {
    mpi_master();
  }
  /*** MPI worker ***/
  else if (m_prank != 0) {
    mpi_worker(output_img);
  }

#elif defined(PARALLEL_MPI) && defined(MPI_SIMPLE)

#ifdef OUTPUT_TIMINGS
  MPI_Barrier(m_communicator);
  auto start = MPI_Wtime();
#endif

  compute_set();

#ifdef OUTPUT_TIMINGS
  auto end = MPI_Wtime();
  if (m_prank == 0) {
    std::cout << m_global_nx << " " << end-start << std::endl;
  }
#endif

  if (output_img) {
    m_pdumper->dump(1, 1);
  }

#else

#ifdef OUTPUT_TIMINGS
  auto start = clk::now();
#endif

  compute_set();

#ifdef OUTPUT_TIMINGS
  auto end = clk::now();
  second compute_time = end - start;
  std::cout << m_global_nx << " " << compute_time.count() << std::endl;
#endif


  if (output_img) {
    m_pdumper->dump(0, 0);
  }
#endif
}


#ifdef PARALLEL_MPI
void Mandelbrot::init_mpi_simple() {

  std::vector<int> locals = get_row_def(m_prank, 
                                        m_global_nx, m_global_ny, m_psize);

  // divide the grid on m_psize rows of height m_local_nx
  m_local_nx = locals[0];
  m_local_ny = locals[1];
  // offsets
  m_local_offset_x = locals[2];
  m_local_offset_y = locals[3];

#ifdef VERBOSE
  std::cerr << m_prank << " " 
            << m_global_nx << " " << m_global_ny << " " 
            << m_local_nx << " " << m_local_ny << " " 
            << m_local_offset_x << " " << m_local_offset_y << std::endl;
#endif

  // resizing the grid
  m_mandel_set.resize(m_local_nx, m_local_ny);
}

void Mandelbrot::mpi_master() {
  int w_prank; // worker prank
  // int row_count;
  int row_idx, n_busy;
  MPI_Status status;
  int w_feeback = 1;

  // number of workers
  int n_workers = m_psize - 1;

  // define a buffer for local sizes
  std::vector<int> buf_locals(4);

  if (n_workers > m_n_rows)
    throw std::invalid_argument("Too much workers processors !");

  std::cerr << "> " << n_workers << " workers" << std::endl;

  /* compute local sizes for each worker + put them in lists to send */
  // initial sendings
  n_busy = 0;
  row_idx = 0;
  for (w_prank = 1; w_prank <= n_workers; w_prank++) {
    // get corresponding local sizes
    buf_locals = get_row_def(row_idx, m_global_nx, m_global_ny, m_n_rows);

    // and send them to worker
    int work_tag = (w_prank == 1 ? FIRST_WORK_TAG : WORK_TAG);
    MPI_Send(&buf_locals[0], 4, MPI_INT, w_prank, work_tag, m_communicator);

    n_busy++;
  }

  int next_row_idx = n_workers;

  /* 'infinite' loop to feed workers with new rows */
  for (;;) {

    MPI_Recv(&w_feeback, 1, MPI_INT, 
             MPI_ANY_SOURCE, FEEBACK_TAG, m_communicator, &status);
    w_prank = status.MPI_SOURCE;

    std::cerr << "--> next row index = " << next_row_idx << std::endl;

    /* send new work */
    if (next_row_idx < m_n_rows) {

      // get corresponding local sizes
      buf_locals = get_row_def(next_row_idx, m_global_nx, m_global_ny, m_n_rows);
      // and send them to worker
      MPI_Send(&buf_locals[0], 4, MPI_INT, w_prank, WORK_TAG, m_communicator);

      if (m_prank == TEST_RANK) std::cerr << "Send buf_locals" << std::endl;

      next_row_idx++;
    }

    /* or tell to quit */
    else {
      int end_tag = (n_busy == 1 ? LAST_END_TAG : END_TAG);
      MPI_Send(&buf_locals[0], 4, MPI_INT, w_prank, end_tag, m_communicator);
      std::cerr << "Send END_TAG to rank " << w_prank << std::endl;
      n_busy--;
    }

    std::cerr << "n_busy " << n_busy << std::endl;

    if (n_busy == 0) {
      std::cerr << "ALL BUSY WORKERS ARE NOTIFIED," << std::endl;
      break; // when no more busy workers
    }

  } /* end infinite loop */
}

void Mandelbrot::mpi_worker(bool output_img) {
  MPI_Status status;
  int w_feeback = 1;

  // define a buffer for local sizes
  std::vector<int> buf_locals(4);

  // set the dumper communicator so that only workers can write the image
  m_pdumper->set_mpi_communicator(m_MW_communicator);

  MPI_Barrier(m_MW_communicator);

  /* infinite loop */
  for (;;) {
    // receive buffer from master containing local problem sizes
    MPI_Recv(&buf_locals[0], 4, MPI_INT, 
             0, MPI_ANY_TAG, m_communicator, &status);

    if (status.MPI_TAG == FIRST_WORK_TAG) {
      std::cerr << m_prank << " first tag !" << std::endl;
      m_pdumper->open_mpi_file();
    }
    else if (status.MPI_TAG == LAST_END_TAG) {
      std::cerr << m_prank << " last tag !" << std::endl;
      m_pdumper->close_mpi_file();
    }

    if ((status.MPI_TAG == LAST_END_TAG) || (status.MPI_TAG == END_TAG)) {
      std::cerr << m_prank << " order to quit ..." << std::endl;
      MPI_Barrier(m_MW_communicator);
      std::cerr << "... byebye rank " << m_prank << std::endl;
      break;
    }

    if (m_prank == TEST_RANK) std::cerr << "Recv buf_locals" << std::endl;

    // unpack variables and update member variables
    m_local_nx       = buf_locals[0];
    m_local_ny       = buf_locals[1];
    m_local_offset_x = buf_locals[2];
    m_local_offset_y = buf_locals[3];

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

    if (output_img) {
      m_pdumper->dump_manual(m_local_offset_x, m_global_nx);
    }

    /* finally send message to tell master that worker is ready for new work */
    MPI_Send(&w_feeback, 1, MPI_INT, 0, FEEBACK_TAG, m_communicator);

    if (m_prank == TEST_RANK) std::cerr << "Send w_feeback" << std::endl;

  } /* end infinite loop */

}

#endif /* PARALLEL_MPI */


void Mandelbrot::compute_set() {
  int ix, iy;

  /* OpenMP parallelization is set here */
#ifdef PARALLEL_OPENMP
  #pragma omp parallel for private(ix, iy) schedule(dynamic)
#endif

  for (ix = 0; ix < m_local_nx; ix++) {

    for (iy = 0; iy < m_local_ny; iy++) {

      compute_pix(ix, iy); /* compute current pixel */

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
  mset(ix, iy) = solve_recursive(cx, cy, z0x, z0y);
}

dfloat Mandelbrot::solve_recursive(dfloat cx, dfloat cy, 
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


std::vector<int> Mandelbrot::get_row_def(int row_idx, int nx, int ny, 
                                         int n_rows) {
  std::vector<int> sizes(4);

  int row_nx = nx / n_rows + (row_idx < nx % n_rows ? 1 : 0);
  int row_ny = ny;
  int row_offset_x = (nx / n_rows) * row_idx + 
                     (row_idx < nx % n_rows ? row_idx : nx % n_rows);
  int row_offset_y = 0;

  sizes[0] = row_nx;
  sizes[1] = row_ny;
  sizes[2] = row_offset_x;
  sizes[3] = row_offset_y;
  return sizes;
}

