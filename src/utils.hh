#include <vector>
#ifdef PARALLEL_MPI
#include <mpi.h>
#endif

// some utility functions

bool all_zeros(std::vector<int> v);

int sum_vector(std::vector<int> v);

void print_vector(std::vector<int> v);
