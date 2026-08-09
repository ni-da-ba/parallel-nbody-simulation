#include "simulation.hpp"

#include <cmath>
#include <iostream>

namespace {
  constexpr double G = 1.0;
  constexpr double SOFTENING = 1e-9;
}

void compute_accelerations(std::vector<Particle>& particles) {
  for (auto& p : particles) {
    p.acceleration = Vec3(0.0, 0.0, 0.0);
  }

  const std::size_t n = particles.size();

  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j < n; ++j) {
      if (i == j) {
	continue;
      }

      Vec3 displacement = particles[j].position - particles[i].position;
      double dist_sq = displacement.norm_squared() + SOFTENING;
      double dist = std::sqrt(dist_sq);
      double inv_dist_cube = 1.0 / (dist_sq * dist);

      double factor = G * particles[j].mass * inv_dist_cube;
      particles[i].acceleration += displacement * factor;
    }
  }
}

void step_particles(std::vector<Particle>& particles, double dt) {
  compute_accelerations(particles);

  for (auto& p : particles) {
    p.velocity += p.acceleration * dt;
    p.position += p.velocity * dt;
  }
}

void compute_accelerations_from_global(
				       std::vector<Particle>& local_particles,
				       const std::vector<Particle>& global_particles
				       ) {
  for (auto& p : local_particles) {
    p.acceleration = Vec3(0.0, 0.0, 0.0);
  }

  for (auto& local_particle : local_particles) {
    for (const auto& global_particle : global_particles) {
      Vec3 displacement = global_particle.position - local_particle.position;
      double dist_sq = displacement.norm_squared() + SOFTENING;

      if (dist_sq <= SOFTENING) {
	continue;
      }

      double dist = std::sqrt(dist_sq);
      double inv_dist_cube = 1.0 / (dist_sq * dist);

      double factor = G * global_particle.mass * inv_dist_cube;
      local_particle.acceleration += displacement * factor;
    }
  }
}

void step_particles_from_global(
				std::vector<Particle>& local_particles,
				const std::vector<Particle>& global_particles,
    double dt
				) {
  compute_accelerations_from_global(local_particles, global_particles);

  for (auto& p : local_particles) {
    p.velocity += p.acceleration * dt;
    p.position += p.velocity * dt;
  }
}

double compute_total_mass(const std::vector<Particle>& particles) {
  double total_mass = 0.0;
  for (const auto& p : particles) {
    total_mass += p.mass;
  }
  return total_mass;
}

Vec3 compute_total_momentum(const std::vector<Particle>& particles) {
  Vec3 total_momentum(0.0, 0.0, 0.0);
  for (const auto& p : particles) {
    total_momentum += p.velocity * p.mass;
  }
  return total_momentum;
}

void print_system_statistics(const std::vector<Particle>& particles, int step) {
  double total_mass = compute_total_mass(particles);
  Vec3 total_momentum = compute_total_momentum(particles);

  std::cout
    << "Stats at step " << step
    << ": total_mass=" << total_mass
    << ", total_momentum=("
    << total_momentum.x << ", "
    << total_momentum.y << ", "
    << total_momentum.z << ")\n";
}

void print_particles(const std::vector<Particle>& particles, int step) {
  std::cout << "Step " << step << '\n';
  for (std::size_t i = 0; i < particles.size(); ++i) {
    const auto& p = particles[i];
    std::cout
      << "Particle " << i
      << " pos=(" << p.position.x << ", " << p.position.y << ", " << p.position.z << ")"
      << " vel=(" << p.velocity.x << ", " << p.velocity.y << ", " << p.velocity.z << ")"
      << '\n';
  }
  print_system_statistics(particles, step);
  std::cout << '\n';
}
