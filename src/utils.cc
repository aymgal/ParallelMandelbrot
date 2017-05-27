#include <algorithm> // std::all_of
#include <iostream>
#include <numeric>   // std::accumulate
#ifdef PARALLEL_MPI
#include <mpi.h>
#endif

#include "utils.hh"
#include "gvars.hh"

// some utility functions

bool all_zeros(std::vector<int> v) {
  return std::all_of(v.begin(), v.end(), [](int i) { return i == 0; });
}

int sum_vector(std::vector<int> v) {
  return std::accumulate(v.begin(), v.end(), 0);
}

void print_vector(std::vector<int> v) {
  for (auto i = v.begin(); i != v.end(); ++i)
    std::cerr << *i << " ";
}