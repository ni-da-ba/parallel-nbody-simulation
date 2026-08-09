#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "distribution.hpp"
#include "octree.hpp"
#include "particle.hpp"
#include "particle_io.hpp"
#include "simulation.hpp"
#include "vec3.hpp"

namespace {

int failures = 0;

void check(bool condition, const std::string& message) {
  if (!condition) {
    std::cerr << "FAIL: " << message << '\n';
    ++failures;
  }
}

bool near(double lhs, double rhs, double tolerance = 1e-12) {
  return std::abs(lhs - rhs) <= tolerance;
}

void check_vec(const Vec3& actual, const Vec3& expected, const std::string& message,
               double tolerance = 1e-12) {
  check(near(actual.x, expected.x, tolerance) &&
            near(actual.y, expected.y, tolerance) &&
            near(actual.z, expected.z, tolerance),
        message);
}

void test_vec3() {
  const Vec3 a(1.0, -2.0, 3.0);
  const Vec3 b(4.0, 5.0, -6.0);
  check_vec(a + b, Vec3(5.0, 3.0, -3.0), "Vec3 addition");
  check_vec(a - b, Vec3(-3.0, -7.0, 9.0), "Vec3 subtraction");
  check_vec(a * 2.0, Vec3(2.0, -4.0, 6.0), "Vec3 scaling");
  check(near(a.norm_squared(), 14.0), "Vec3 squared norm");
}

void test_distribution() {
  check(all_local_counts(10, 3) == std::vector<int>({4, 3, 3}),
        "MPI particle counts include the remainder");
  check(all_local_offsets(10, 3) == std::vector<int>({0, 4, 7}),
        "MPI particle offsets are contiguous");
  check(all_local_counts(2, 4) == std::vector<int>({1, 1, 0, 0}),
        "More ranks than particles are handled");
}

void test_particle_serialization() {
  const std::vector<Particle> particles = {
      Particle(Vec3(1.0, 2.0, 3.0), Vec3(4.0, 5.0, 6.0), 7.0),
      Particle(Vec3(-1.0, -2.0, -3.0), Vec3(-4.0, -5.0, -6.0), 0.5),
  };
  const auto unpacked = unpack_particles(pack_particles(particles));
  check(unpacked.size() == particles.size(), "Particle round-trip count");
  check_vec(unpacked[0].position, particles[0].position,
            "Particle position round trip");
  check_vec(unpacked[1].velocity, particles[1].velocity,
            "Particle velocity round trip");
  check(near(unpacked[1].mass, particles[1].mass),
        "Particle mass round trip");

  bool rejected = false;
  try {
    static_cast<void>(unpack_particles({1.0, 2.0}));
  } catch (const std::runtime_error&) {
    rejected = true;
  }
  check(rejected, "Malformed particle buffers are rejected");
}

void test_dynamics() {
  std::vector<Particle> particles = {
      Particle(Vec3(-1.0, 0.0, 0.0), Vec3(0.0, 1.0, 0.0), 2.0),
      Particle(Vec3(1.0, 0.0, 0.0), Vec3(0.0, -2.0, 0.0), 1.0),
  };

  check(near(compute_total_mass(particles), 3.0), "Total mass");
  check_vec(compute_total_momentum(particles), Vec3(0.0, 0.0, 0.0),
            "Total momentum");

  compute_accelerations(particles);
  check(particles[0].acceleration.x > 0.0,
        "Left particle accelerates toward the right particle");
  check(particles[1].acceleration.x < 0.0,
        "Right particle accelerates toward the left particle");
  check(near(2.0 * particles[0].acceleration.x +
                 particles[1].acceleration.x,
             0.0),
        "Direct forces preserve equal-and-opposite momentum change");
}

void test_octree() {
  const std::vector<Particle> particles = {
      Particle(Vec3(-1.0, 0.0, 0.0), Vec3(), 1.0),
      Particle(Vec3(1.0, 0.0, 0.0), Vec3(), 2.0),
      Particle(Vec3(0.0, 1.0, 0.5), Vec3(), 0.5),
  };

  const Cube bounds = make_bounding_cube(particles);
  for (const auto& particle : particles) {
    check(bounds.contains(particle.position),
          "Bounding cube contains every particle");
  }

  OctreeNode root(bounds);
  for (int i = 0; i < static_cast<int>(particles.size()); ++i) {
    root.insert(i, particles);
  }
  root.compute_mass_properties(particles);
  check(near(root.total_mass(), 3.5), "Octree aggregates total mass");

  auto direct = particles;
  compute_accelerations(direct);
  const Vec3 tree_acceleration = root.compute_acceleration_on_particle(
      0, particles, 1e-12, 1e-9, 1.0);
  check_vec(tree_acceleration, direct[0].acceleration,
            "Exact octree traversal agrees with direct acceleration", 1e-11);

  bool rejected_empty = false;
  try {
    static_cast<void>(make_bounding_cube({}));
  } catch (const std::invalid_argument&) {
    rejected_empty = true;
  }
  check(rejected_empty, "Empty octree inputs are rejected explicitly");

  const std::vector<Particle> coincident = {
      Particle(Vec3(1.0, 1.0, 1.0), Vec3(), 1.0),
      Particle(Vec3(1.0, 1.0, 1.0), Vec3(), 2.0),
  };
  OctreeNode coincident_root(Cube(Vec3(1.0, 1.0, 1.0), 1.0));
  coincident_root.insert(0, coincident);
  bool rejected_coincident = false;
  try {
    coincident_root.insert(1, coincident);
  } catch (const std::runtime_error&) {
    rejected_coincident = true;
  }
  check(rejected_coincident,
        "Coincident particles fail explicitly instead of recursing indefinitely");
}

}  // namespace

int main() {
  test_vec3();
  test_distribution();
  test_particle_serialization();
  test_dynamics();
  test_octree();

  if (failures != 0) {
    std::cerr << failures << " test(s) failed\n";
    return 1;
  }

  std::cout << "All core tests passed\n";
  return 0;
}
