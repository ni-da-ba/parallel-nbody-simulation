#ifndef PARTICLE_IO_HPP
#define PARTICLE_IO_HPP

#include <vector>

#include "particle.hpp"

constexpr int PARTICLE_DOUBLES = 7;

std::vector<double> pack_particles(const std::vector<Particle>& particles);
std::vector<Particle> unpack_particles(const std::vector<double>& buffer);

#endif
