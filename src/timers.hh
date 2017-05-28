#ifdef PARALLEL_MPI
#include <mpi.h>
#endif

#include "gvars.hh"


class Timer {
public:
  Timer() {}

  virtual void start_chrono() = 0;
  virtual void stop_chrono()   = 0;
  virtual double get_timing() = 0;

protected:
  double m_timing;
};


class TimerSTD : public Timer {
public:
  TimerSTD();

  virtual void start_chrono();
  virtual void stop_chrono();
  virtual double get_timing();

private:
  tp_type m_start_time;
  tp_type m_end_time;
};


#ifdef PARALLEL_MPI

class TimerMPI : public Timer {
public:
  TimerMPI(MPI_Comm comm);

  virtual void start_chrono();
  virtual void stop_chrono();
  virtual double get_timing();

private:
  double m_start_time;
  double m_end_time;
  MPI_Comm m_communicator;
};

#endif


