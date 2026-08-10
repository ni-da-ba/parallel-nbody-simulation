#ifndef OCTREE_HPP
#define OCTREE_HPP

#include <array>
#include <memory>
#include <vector>

#include "particle.hpp"
#include "vec3.hpp"

struct Cube {
  Vec3 center;
  double half_size;

  Cube();
  Cube(const Vec3& center_, double half_size_);

  bool contains(const Vec3& point) const;
  Cube child_cube(int octant) const;
};

class OctreeNode {
 public:
  explicit OctreeNode(const Cube& region_);

  bool is_leaf() const;
  void insert(int particle_index, const std::vector<Particle>& particles);
  void compute_mass_properties(const std::vector<Particle>& particles);
  void print(const std::vector<Particle>& particles, int depth = 0) const;

  Vec3 compute_acceleration_on_particle(int target_particle_index,
                                        const std::vector<Particle>& particles, double theta,
                                        double softening, double G) const;

  const Cube& region() const;
  double total_mass() const;
  const Vec3& center_of_mass() const;
  int particle_index() const;

 private:
  Cube region_;
  double total_mass_;
  Vec3 center_of_mass_;
  int particle_index_;
  std::array<std::unique_ptr<OctreeNode>, 8> children_;

  int get_octant(const Vec3& point) const;
  void subdivide();
};

Cube make_bounding_cube(const std::vector<Particle>& particles);

#endif
