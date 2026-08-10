#ifndef SIMULATION_HPP
#define SIMULATION_HPP

#include <vector>

#include "particle.hpp"

void compute_accelerations(std::vector<Particle>& particles);
void step_particles(std::vector<Particle>& particles, double dt);

void compute_accelerations_from_global(std::vector<Particle>& local_particles,
                                       const std::vector<Particle>& global_particles);

void step_particles_from_global(std::vector<Particle>& local_particles,
                                const std::vector<Particle>& global_particles, double dt);

double compute_total_mass(const std::vector<Particle>& particles);
Vec3 compute_total_momentum(const std::vector<Particle>& particles);
void print_system_statistics(const std::vector<Particle>& particles, int step);
void print_particles(const std::vector<Particle>& particles, int step);

#endif
