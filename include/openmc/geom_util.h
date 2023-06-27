#ifndef OPENMC_GEOM_UTIL_H
#define OPENMC_GEOM_UTIL_H

#include <array>

#include "openmc/position.h"

namespace openmc {

bool plucker_ray_tri_intersect(const std::array<Position, 3> vertices,
  const Position& origin, const Position& direction, double& dist_out,
  const double* neg_ray_len = nullptr, const int* orientation = nullptr);
}

#endif