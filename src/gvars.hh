// global variables

#ifndef GVARS_HH
#define GVARS_HH

#include <chrono>

#ifdef DOUBLE_PRECISION
typedef double dfloat;
#else
typedef float dfloat;
#endif

typedef std::chrono::high_resolution_clock::time_point tp_type;
typedef std::chrono::high_resolution_clock clk;
typedef std::chrono::duration<double> seconds;

#endif /* GVARS_HH */