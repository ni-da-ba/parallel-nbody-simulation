#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include "vec3.hpp"

struct Particle {
  Vec3 position;
  Vec3 velocity;
  Vec3 acceleration;
  double mass;

  Particle() : position(), velocity(), acceleration(), mass(1.0) {}

  Particle(const Vec3& pos, const Vec3& vel, double m)
      : position(pos), velocity(vel), acceleration(), mass(m) {}
};

#endif
