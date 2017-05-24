// global variables

#ifndef GVARS_HH
#define GVARS_HH

#include <chrono>

typedef std::chrono::high_resolution_clock clk;
typedef std::chrono::duration<double> second;

// define limits of the image
#define XMIN -2.0
#define XMAX  1.0
#define YMIN -1.5
#define YMAX  1.5

#ifdef DOUBLE_PRECISION
typedef double dfloat;
#else
typedef float dfloat;
#endif

#endif /* GVARS_HH */