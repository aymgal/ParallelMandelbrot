#ifdef PARALLEL_MPI
#include <mpi.h>
#endif

#include "gvars.hh"

class TimerSTD {
public:
  TimerSTD()
    : m_start_time(clk::now()), m_end_time(clk::now()) {}

  void start_chrono() {
    m_start_time = clk::now(); 
  }

  void end_chrono() {
    m_end_time = clk::now(); 
  }

  seconds get_timing() {
    m_timing = m_end_time - m_start_time;
    return m_timing;
  }

private:
  tp_type m_start_time;
  tp_type m_end_time;
  seconds m_timing;
};


#ifdef PARALLEL_MPI

class TimerMPI {
public:
  TimerMPI(MPI_Comm comm) 
    : m_communicator(comm), m_start_time(0.0), m_end_time(0.0) {}

  void start_chrono() {
    MPI_Barrier(m_communicator);
    m_start_time = MPI_Wtime(); 
  }

  void end_chrono() {
    MPI_Barrier(m_communicator);
    m_end_time = MPI_Wtime(); 
  }

  double get_timing() {
    m_timing = m_end_time - m_start_time;
    return m_timing;
  }

private:
  MPI_Comm m_communicator;
  double m_start_time;
  double m_end_time;
  double m_timing;
};

#endif


