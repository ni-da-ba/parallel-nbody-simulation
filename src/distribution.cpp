#include "distribution.hpp"

int local_particle_count(int global_n, int rank, int size) {
  const int base = global_n / size;
  const int remainder = global_n % size;
  return base + (rank < remainder ? 1 : 0);
}

int local_particle_offset(int global_n, int rank, int size) {
  const int base = global_n / size;
  const int remainder = global_n % size;

  if (rank < remainder) {
    return rank * (base + 1);
  }

  return remainder * (base + 1) + (rank - remainder) * base;
}

std::vector<int> all_local_counts(int global_n, int size) {
  std::vector<int> counts(size);
  for (int rank = 0; rank < size; ++rank) {
    counts[rank] = local_particle_count(global_n, rank, size);
  }
  return counts;
}

std::vector<int> all_local_offsets(int global_n, int size) {
  std::vector<int> offsets(size);
  for (int rank = 0; rank < size; ++rank) {
    offsets[rank] = local_particle_offset(global_n, rank, size);
  }
  return offsets;
}
