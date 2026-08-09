#include "particle_io.hpp"

#include <stdexcept>

std::vector<double> pack_particles(const std::vector<Particle>& particles) {
  std::vector<double> buffer;
  buffer.reserve(static_cast<std::size_t>(PARTICLE_DOUBLES) * particles.size());

  for (const auto& p : particles) {
    buffer.push_back(p.position.x);
    buffer.push_back(p.position.y);
    buffer.push_back(p.position.z);

    buffer.push_back(p.velocity.x);
    buffer.push_back(p.velocity.y);
    buffer.push_back(p.velocity.z);

    buffer.push_back(p.mass);
  }

  return buffer;
}

std::vector<Particle> unpack_particles(const std::vector<double>& buffer) {
  if (buffer.size() % PARTICLE_DOUBLES != 0) {
    throw std::runtime_error("Particle buffer size is not divisible by PARTICLE_DOUBLES.");
  }

  std::vector<Particle> particles;
  particles.reserve(buffer.size() / PARTICLE_DOUBLES);

  for (std::size_t i = 0; i < buffer.size(); i += PARTICLE_DOUBLES) {
    Vec3 position(buffer[i], buffer[i + 1], buffer[i + 2]);
    Vec3 velocity(buffer[i + 3], buffer[i + 4], buffer[i + 5]);
    double mass = buffer[i + 6];

    particles.emplace_back(position, velocity, mass);
  }

  return particles;
}
