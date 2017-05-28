#ifdef PARALLEL_MPI
#include <mpi.h>
#endif

#ifndef DUMPERS_HH
#define DUMPERS_HH

class Grid;

class Dumper {
public:
  explicit Dumper(const Grid & grid)
      : m_grid(grid), m_min(0.0), m_max(255.0) {}

  virtual void dump(int arg1, int arg2) = 0;

  void set_min(float min);
  void set_max(float max);

  float min() const;
  float max() const;

protected:
  const Grid & m_grid;
  float m_min, m_max;
};

class DumperASCII : public Dumper {
public:
  explicit DumperASCII(const Grid & grid)
      : Dumper(grid) {}

  virtual void dump(int arg1, int arg2);
};

#ifdef PARALLEL_MPI
class DumperBinary : public Dumper {
public:
  explicit DumperBinary(const Grid & grid, MPI_Comm comm)
      : Dumper(grid), m_communicator(comm), m_header(false) {}

  virtual void dump(int arg1, int arg2);

  void dump_manual(int arg1, int arg2, int offset_h, int total_h);

  void set_mpi_communicator(MPI_Comm new_comm);

private:
  MPI_Comm m_communicator;
  bool m_header;
};
#endif

#endif /* DUMPERS_HH */
