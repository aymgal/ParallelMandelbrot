#include <sstream>
#ifdef PARALLEL_MPI
#include <mpi.h>
#endif

#ifndef DUMPERS_HH
#define DUMPERS_HH

class Grid;

class Dumper {
public:
#ifdef PARALLEL_MPI
  explicit Dumper(const Grid & grid, MPI_Comm comm)
      : m_grid(grid), m_min(0.0), m_max(255.0), m_communicator(comm) {}
#else
  explicit Dumper(const Grid & grid)
      : m_grid(grid), m_min(0.0), m_max(255.0) {}
#endif

  virtual void dump(int arg1, int arg2) = 0;

  void set_min(float min);
  void set_max(float max);

  float min() const;
  float max() const;

protected:
  const Grid & m_grid;
  float m_min, m_max;
  std::stringstream m_sfilename;

#ifdef PARALLEL_MPI
  MPI_Comm m_communicator;
#endif

};

class DumperASCII : public Dumper {
public:
#ifdef PARALLEL_MPI
  explicit DumperASCII(const Grid & grid, MPI_Comm comm)
      : Dumper(grid, comm) {
    m_sfilename << "out_default.pgm";
  }
#else
  explicit DumperASCII(const Grid & grid)
      : Dumper(grid) {}
#endif

  virtual void dump(int arg1, int arg2);
};

#ifdef PARALLEL_MPI
class DumperBinary : public Dumper {
public:
  explicit DumperBinary(const Grid & grid, MPI_Comm comm)
      : Dumper(grid, comm), m_header_written(false) {
    m_sfilename << "out_default.bmp";
  }

  virtual void dump(int arg1, int arg2);

  void dump_manual(int offset_h, int total_h);

  void set_mpi_communicator(MPI_Comm new_comm);
  void open_mpi_file();
  void close_mpi_file();

private:
  bool m_header_written;
  MPI_File m_fh;
};
#endif

#endif /* DUMPERS_HH */
