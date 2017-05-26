#include <array>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "dumpers.hh"
#include "grid.hh"
#include "gvars.hh"

#define COLOR_LEVELS 100.0


void Dumper::set_min(float min) {
  m_min = min;
}

void Dumper::set_max(float max) {
  m_max = max;
}

float Dumper::min() const {
  return m_min;
}

float Dumper::max() const {
  return m_max;
}

void DumperASCII::dump(int arg1, int arg2) {
  std::ofstream fout;

  int m = m_grid.nx();
  int n = m_grid.ny();

  set_min(0.0);
  // set_max(m_grid.get_max());
  set_max(COLOR_LEVELS);

  m_sfilename << "out_" << arg1 << "_" << arg2 << ".pgm";

  fout.open(m_sfilename.str());

  fout << "P2" << std::endl 
       << "# CREATOR: Mandelbrot serial program" << std::endl;
  fout << m << " " << n << std::endl;
  fout << 255 << std::endl;

  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      int v = 255. * (m_grid(j, i) - m_min) / (m_max - m_min);
      v = std::min(v, 255);
      fout << v << std::endl;
    }
  }

  fout.close();
}


#ifdef PARALLEL_MPI

void DumperBinary::dump(int arg1, int arg2) {

  int prank, psize;
  MPI_Comm_rank(m_communicator, &prank);
  MPI_Comm_size(m_communicator, &psize);

  m_sfilename << "out_" << arg1 << "_" << arg2 << ".bmp";

  int h = m_grid.nx();
  int w = m_grid.ny();

  set_min(0.0);
#if defined(PARALLEL_MPI) && defined(COLOR_PRANK)
  set_max(psize);
#else
  // set_max(m_grid.get_max());
  set_max(COLOR_LEVELS);
#endif

  // Gathering the size of every processors, this could be done as in the
  // constructor of the Simulation instead
  std::vector<int> size_per_proc(psize);
  MPI_Allgather(&h, 1, MPI_INT, size_per_proc.data(), 1, MPI_INT,
                m_communicator);

  // determining the local offset
  int offset_h = 0;
  for (int i = 0; i < prank; ++i) {
    offset_h += size_per_proc[i];
  }

  int total_h = offset_h;
  for (int i = prank; i < psize; ++i) {
    total_h += size_per_proc[i];
  }

  int row_size = 3 * w;
  // if the file width (3*w) is not a multiple of 4 adds enough bytes to make it
  // a multiple of 4
  int padding = (4 - (row_size) % 4) % 4;
  row_size += padding;

  int filesize = 54 + row_size * total_h;

  std::vector<char> img(row_size * h);
  std::fill(img.begin(), img.end(), 0);

  for (int i = 0; i < h; i++) {
    for (int j = 0; j < w; j++) {
      float v = ((m_grid(i, j) - m_min) / (m_max - m_min));

      float r = v * 255.0; // R channel
      float g = v * 255.0; // G channel
      float b = v * 255.0; // B channel

      r = std::min(r, 255.f);
      g = std::min(g, 255.f);
      b = std::min(b, 255.f);

      img[row_size * i + 3 * j + 2] = r;
      img[row_size * i + 3 * j + 1] = g;
      img[row_size * i + 3 * j + 0] = b;
    }
  }

  std::array<char, 14> bmpfileheader = {{'B', 'M', 0, 0,  0, 0, 0,
                                        0,   0,   0, 54, 0, 0, 0}};
  std::array<char, 40> bmpinfoheader = {{40, 0, 0, 0, 0, 0, 0,  0,
                                        0,  0, 0, 0, 1, 0, 24, 0}};

  bmpfileheader[2] = filesize;
  bmpfileheader[3] = filesize >> 8;
  bmpfileheader[4] = filesize >> 16;
  bmpfileheader[5] = filesize >> 24;

  bmpinfoheader[4] = w;
  bmpinfoheader[5] = w >> 8;
  bmpinfoheader[6] = w >> 16;
  bmpinfoheader[7] = w >> 24;
  bmpinfoheader[8] = total_h;
  bmpinfoheader[9] = total_h >> 8;
  bmpinfoheader[10] = total_h >> 16;
  bmpinfoheader[11] = total_h >> 24;
  bmpinfoheader[20] = (filesize - 54);
  bmpinfoheader[21] = (filesize - 54) >> 8;
  bmpinfoheader[22] = (filesize - 54) >> 16;
  bmpinfoheader[23] = (filesize - 54) >> 24;

  MPI_Status status;

  open_mpi_file();

  // defining the size of the file
  MPI_File_set_size(m_fh, filesize);

  // write the header
  if (!m_header_written) {
    MPI_File_write_at(m_fh, 0, bmpfileheader.data(), 14, MPI_CHAR, &status);
    MPI_File_write_at(m_fh, 14, bmpinfoheader.data(), 40, MPI_CHAR, &status);
    m_header_written = true;
  }

  int offset = 54 + row_size * offset_h;

  // We also could write that data with a write_at, the view is just to show
  // different possibilities
  MPI_File_set_view(m_fh, offset, MPI_CHAR, MPI_CHAR, "native", MPI_INFO_NULL);

  MPI_File_write(m_fh, img.data(), img.size(), MPI_CHAR, &status);
  close_mpi_file();
}


void DumperBinary::dump_manual(int offset_h, int total_h) {

  int prank, psize;
  MPI_Comm_rank(m_communicator, &prank);
  MPI_Comm_size(m_communicator, &psize);

  // m_sfilename << "out_" << "manual" << ".bmp";

  int h = m_grid.nx();
  int w = m_grid.ny();

  set_min(0.0);
#if defined(PARALLEL_MPI) && defined(COLOR_PRANK)
  set_max(psize);
#else
  // set_max(m_grid.get_max());
  set_max(COLOR_LEVELS);
#endif

  int row_size = 3 * w;
  // if the file width (3*w) is not a multiple of 4 adds enough bytes to make it
  // a multiple of 4
  int padding = (4 - (row_size) % 4) % 4;
  row_size += padding;

  int filesize = 54 + row_size * total_h;

  std::vector<char> img(row_size * h);
  std::fill(img.begin(), img.end(), 0);

  for (int i = 0; i < h; i++) {
    for (int j = 0; j < w; j++) {
      float v = ((m_grid(i, j) - m_min) / (m_max - m_min));

      float r = v * 255.0; // R channel
      float g = v * 255.0; // G channel
      float b = v * 255.0; // B channel

      r = std::min(r, 255.f);
      g = std::min(g, 255.f);
      b = std::min(b, 255.f);

      img[row_size * i + 3 * j + 2] = r;
      img[row_size * i + 3 * j + 1] = g;
      img[row_size * i + 3 * j + 0] = b;
    }
  }

  std::array<char, 14> bmpfileheader = {{'B', 'M', 0, 0,  0, 0, 0,
                                        0,   0,   0, 54, 0, 0, 0}};
  std::array<char, 40> bmpinfoheader = {{40, 0, 0, 0, 0, 0, 0,  0,
                                        0,  0, 0, 0, 1, 0, 24, 0}};

  bmpfileheader[2] = filesize;
  bmpfileheader[3] = filesize >> 8;
  bmpfileheader[4] = filesize >> 16;
  bmpfileheader[5] = filesize >> 24;

  bmpinfoheader[4] = w;
  bmpinfoheader[5] = w >> 8;
  bmpinfoheader[6] = w >> 16;
  bmpinfoheader[7] = w >> 24;
  bmpinfoheader[8] = total_h;
  bmpinfoheader[9] = total_h >> 8;
  bmpinfoheader[10] = total_h >> 16;
  bmpinfoheader[11] = total_h >> 24;
  bmpinfoheader[20] = (filesize - 54);
  bmpinfoheader[21] = (filesize - 54) >> 8;
  bmpinfoheader[22] = (filesize - 54) >> 16;
  bmpinfoheader[23] = (filesize - 54) >> 24;

  MPI_Status status;

  // std::cerr << prank+1 << " FUCK-1" << std::endl;
  // MPI_File_open(m_communicator, m_sfilename.str().c_str(),
  //               MPI_MODE_WRONLY | MPI_MODE_CREATE,
  //               MPI_INFO_NULL, &m_fh);
  // std::cerr << prank+1 << " FUCK-2" << std::endl;

  // defining the size of the file
  MPI_File_set_size(m_fh, filesize);

  // rank 0 writes the header
  // if (prank == 0) {
  if (!m_header_written) {
    MPI_File_write_at(m_fh, 0, bmpfileheader.data(), 14, MPI_CHAR, &status);
    MPI_File_write_at(m_fh, 14, bmpinfoheader.data(), 40, MPI_CHAR, &status);
    m_header_written = true;
  }

  int offset = 54 + row_size * offset_h;

  // We also could write that data with a write_at, the view is just to show
  // different possibilities
  MPI_File_set_view(m_fh, offset, MPI_CHAR, MPI_CHAR, "native", MPI_INFO_NULL);

  MPI_File_write(m_fh, img.data(), img.size(), MPI_CHAR, &status);
  // MPI_File_close(&m_fh);
}

void DumperBinary::set_mpi_communicator(MPI_Comm new_comm) {
  m_communicator = new_comm;
}

void DumperBinary::open_mpi_file() {
  std::cerr << "FUCK 1" << std::endl;
  // opening a file in write and create mode
  MPI_File_open(m_communicator, m_sfilename.str().c_str(),
                MPI_MODE_WRONLY | MPI_MODE_CREATE,
                MPI_INFO_NULL, &m_fh);
  std::cerr << "FUCK 2" << std::endl;
}

void DumperBinary::close_mpi_file() {
  MPI_File_close(&m_fh);
}

#endif
