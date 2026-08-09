#ifndef VEC3_HPP
#define VEC3_HPP

#include <cmath>

struct Vec3 {
  double x;
  double y;
  double z;

  Vec3() : x(0.0), y(0.0), z(0.0) {}
  Vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}

  Vec3 operator+(const Vec3& other) const {
    return Vec3(x + other.x, y + other.y, z + other.z);
  }

  Vec3 operator-(const Vec3& other) const {
    return Vec3(x - other.x, y - other.y, z - other.z);
  }

  Vec3 operator*(double scalar) const {
    return Vec3(x * scalar, y * scalar, z * scalar);
  }

  Vec3 operator/(double scalar) const {
    return Vec3(x / scalar, y / scalar, z / scalar);
  }

  Vec3& operator+=(const Vec3& other) {
    x += other.x;
    y += other.y;
    z += other.z;
    return *this;
  }

  Vec3& operator-=(const Vec3& other) {
    x -= other.x;
    y -= other.y;
    z -= other.z;
    return *this;
  }

  Vec3& operator*=(double scalar) {
    x *= scalar;
    y *= scalar;
    z *= scalar;
    return *this;
  }

  double norm_squared() const {
    return x * x + y * y + z * z;
  }

  double norm() const {
    return std::sqrt(norm_squared());
  }
};

#endif
