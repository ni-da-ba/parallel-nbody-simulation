#ifndef DISTRIBUTION_HPP
#define DISTRIBUTION_HPP

#include <vector>

int local_particle_count(int global_n, int rank, int size);
int local_particle_offset(int global_n, int rank, int size);

std::vector<int> all_local_counts(int global_n, int size);
std::vector<int> all_local_offsets(int global_n, int size);

#endif
