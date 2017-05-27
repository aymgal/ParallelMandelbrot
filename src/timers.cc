#include "timers.hh"

TimerSTD::TimerSTD()
  : Timer(), m_start_time(clk::now()), m_end_time(clk::now()) {}

void TimerSTD::start_chrono() {
  m_start_time = clk::now();
}

void TimerSTD::end_chrono() {
  m_end_time = clk::now();
}

double TimerSTD::get_timing() {
  seconds timing_tmp = m_end_time - m_start_time;
  m_timing = timing_tmp.count();
  return m_timing;
}

#ifdef PARALLEL_MPI

TimerMPI::TimerMPI(MPI_Comm comm) 
    : Timer(), m_start_time(0.0), m_end_time(0.0), m_communicator(comm) {}

void TimerMPI::start_chrono() {
  MPI_Barrier(m_communicator);
  m_start_time = MPI_Wtime();
}

void TimerMPI::end_chrono() {
  MPI_Barrier(m_communicator);
  m_end_time = MPI_Wtime();
}

double TimerMPI::get_timing() {
  m_timing = m_end_time - m_start_time;
  return m_timing;
}

#endif
