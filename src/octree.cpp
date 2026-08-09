#include "octree.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

Cube::Cube() : center(), half_size(0.0) {}

Cube::Cube(const Vec3& center_, double half_size_)
  : center(center_), half_size(half_size_) {}

bool Cube::contains(const Vec3& point) const {
    return
        point.x >= center.x - half_size && point.x <= center.x + half_size &&
      point.y >= center.y - half_size && point.y <= center.y + half_size &&
      point.z >= center.z - half_size && point.z <= center.z + half_size;
}

Cube Cube::child_cube(int octant) const {
  const double offset = half_size / 2.0;

  const double child_x = center.x + ((octant & 1) ? offset : -offset);
  const double child_y = center.y + ((octant & 2) ? offset : -offset);
  const double child_z = center.z + ((octant & 4) ? offset : -offset);

  return Cube(Vec3(child_x, child_y, child_z), offset);
}

OctreeNode::OctreeNode(const Cube& region_)
  : region_(region_),
    total_mass_(0.0),
    center_of_mass_(0.0, 0.0, 0.0),
    particle_index_(-1),
    children_() {}

bool OctreeNode::is_leaf() const {
  for (const auto& child : children_) {
    if (child) {
      return false;
    }
  }
  return true;
}

const Cube& OctreeNode::region() const {
  return region_;
}

double OctreeNode::total_mass() const {
  return total_mass_;
}

const Vec3& OctreeNode::center_of_mass() const {
  return center_of_mass_;
}

int OctreeNode::particle_index() const {
  return particle_index_;
}

int OctreeNode::get_octant(const Vec3& point) const {
  int octant = 0;

  if (point.x >= region_.center.x) {
    octant |= 1;
  }
  if (point.y >= region_.center.y) {
    octant |= 2;
  }
  if (point.z >= region_.center.z) {
    octant |= 4;
  }

  return octant;
}

void OctreeNode::subdivide() {
  for (int octant = 0; octant < 8; ++octant) {
    children_[octant] = std::make_unique<OctreeNode>(region_.child_cube(octant));
  }
}

void OctreeNode::insert(int new_particle_index, const std::vector<Particle>& particles) {
  const Vec3& pos = particles[new_particle_index].position;

  if (!region_.contains(pos)) {
    throw std::runtime_error("Tried to insert particle outside octree region.");
  }

  if (is_leaf() && particle_index_ == -1) {
    particle_index_ = new_particle_index;
    return;
  }

  if (is_leaf()) {
    const int existing_particle_index = particle_index_;

    const Vec3& existing_pos = particles[existing_particle_index].position;
    if (existing_pos.x == pos.x && existing_pos.y == pos.y &&
        existing_pos.z == pos.z) {
      throw std::runtime_error(
          "Barnes-Hut tree cannot separate coincident particle positions.");
    }

    particle_index_ = -1;
    subdivide();

    const int existing_octant = get_octant(particles[existing_particle_index].position);
    children_[existing_octant]->insert(existing_particle_index, particles);
  }

  const int new_octant = get_octant(pos);
  children_[new_octant]->insert(new_particle_index, particles);
}

void OctreeNode::compute_mass_properties(const std::vector<Particle>& particles) {
  total_mass_ = 0.0;
  center_of_mass_ = Vec3(0.0, 0.0, 0.0);

  if (is_leaf()) {
    if (particle_index_ >= 0) {
      total_mass_ = particles[particle_index_].mass;
      center_of_mass_ = particles[particle_index_].position;
    }
    return;
  }

  Vec3 weighted_position_sum(0.0, 0.0, 0.0);

  for (auto& child : children_) {
    if (!child) {
      continue;
    }

    child->compute_mass_properties(particles);

    total_mass_ += child->total_mass();
    weighted_position_sum += child->center_of_mass() * child->total_mass();
  }

  if (total_mass_ > 0.0) {
    center_of_mass_ = weighted_position_sum / total_mass_;
  }
}

Vec3 OctreeNode::compute_acceleration_on_particle(
						  int target_particle_index,
						  const std::vector<Particle>& particles,
						  double theta,
						  double softening,
    double G
						  ) const {
  if (total_mass_ <= 0.0) {
    return Vec3(0.0, 0.0, 0.0);
  }

  const Particle& target = particles[target_particle_index];

  if (is_leaf()) {
    if (particle_index_ < 0 || particle_index_ == target_particle_index) {
      return Vec3(0.0, 0.0, 0.0);
    }

    Vec3 displacement = particles[particle_index_].position - target.position;
    double dist_sq = displacement.norm_squared() + softening;
    double dist = std::sqrt(dist_sq);
    double inv_dist_cube = 1.0 / (dist_sq * dist);

    return displacement * (G * particles[particle_index_].mass * inv_dist_cube);
  }

  Vec3 displacement = center_of_mass_ - target.position;
  double dist = displacement.norm();

  if (dist == 0.0) {
    return Vec3(0.0, 0.0, 0.0);
  }

  const double s = 2.0 * region_.half_size;

  if ((s / dist) < theta) {
    double dist_sq = displacement.norm_squared() + softening;
    double softened_dist = std::sqrt(dist_sq);
    double inv_dist_cube = 1.0 / (dist_sq * softened_dist);

    return displacement * (G * total_mass_ * inv_dist_cube);
  }

  Vec3 total_acceleration(0.0, 0.0, 0.0);

  for (const auto& child : children_) {
    if (child) {
      total_acceleration += child->compute_acceleration_on_particle(
								    target_particle_index, particles, theta, softening, G
								    );
    }
  }

  return total_acceleration;
}

void OctreeNode::print(const std::vector<Particle>& particles, int depth) const {
  const std::string indent(static_cast<std::size_t>(depth) * 2, ' ');

  std::cout
    << indent
    << "Node(center=("
    << region_.center.x << ", "
    << region_.center.y << ", "
    << region_.center.z << "), half_size="
    << region_.half_size
    << ", total_mass=" << total_mass_
    << ", center_of_mass=("
    << center_of_mass_.x << ", "
    << center_of_mass_.y << ", "
    << center_of_mass_.z << ")";

  if (is_leaf()) {
    std::cout << ", leaf_particle=" << particle_index_;
    if (particle_index_ >= 0) {
      const auto& p = particles[particle_index_];
      std::cout << " @ pos=("
		<< p.position.x << ", "
		<< p.position.y << ", "
		<< p.position.z << ")";
    }
  }

  std::cout << ")\n";

  if (!is_leaf()) {
    for (const auto& child : children_) {
      if (child) {
	child->print(particles, depth + 1);
      }
    }
  }
}

Cube make_bounding_cube(const std::vector<Particle>& particles) {
  if (particles.empty()) {
    throw std::invalid_argument(
        "Cannot construct an octree bounding cube for zero particles.");
  }

  double min_x = particles[0].position.x;
  double max_x = particles[0].position.x;
  double min_y = particles[0].position.y;
  double max_y = particles[0].position.y;
  double min_z = particles[0].position.z;
  double max_z = particles[0].position.z;

  for (const auto& p : particles) {
    min_x = std::min(min_x, p.position.x);
    max_x = std::max(max_x, p.position.x);
    min_y = std::min(min_y, p.position.y);
    max_y = std::max(max_y, p.position.y);
    min_z = std::min(min_z, p.position.z);
    max_z = std::max(max_z, p.position.z);
  }

  const Vec3 center(
		    0.5 * (min_x + max_x),
		    0.5 * (min_y + max_y),
		    0.5 * (min_z + max_z)
		    );

  const double span_x = max_x - min_x;
  const double span_y = max_y - min_y;
  const double span_z = max_z - min_z;
  double span = std::max(span_x, std::max(span_y, span_z));

  if (span == 0.0) {
    span = 1.0;
  }

  const double half_size = 0.5 * span + 1e-6;
  return Cube(center, half_size);
}
