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

#if PARALLEL_MPI
  virtual void set_communicator(MPI_Comm new_comm) = 0;
#endif

protected:
  const Grid & m_grid;
  float m_min, m_max;
#ifdef PARALLEL_MPI
  MPI_Comm m_communicator;
#endif

};

class DumperASCII : public Dumper {
public:
#ifdef PARALLEL_MPI
  explicit DumperASCII(const Grid & grid, MPI_Comm comm)
      : Dumper(grid, comm) {}
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
      : Dumper(grid, comm) {}

  virtual void dump(int arg1, int arg2);

  virtual void set_communicator(MPI_Comm new_comm);
};
#endif

#endif /* DUMPERS_HH */
