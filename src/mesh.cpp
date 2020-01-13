#include "openmc/mesh.h"

#include <algorithm> // for copy, equal, min, min_element
#include <cstddef> // for size_t
#include <cmath>  // for ceil
#include <memory> // for allocator
#include <string>
#include <gsl/gsl>

#ifdef OPENMC_MPI
#include "mpi.h"
#endif
#include "xtensor/xbuilder.hpp"
#include "xtensor/xeval.hpp"
#include "xtensor/xmath.hpp"
#include "xtensor/xsort.hpp"
#include "xtensor/xtensor.hpp"

#include "openmc/capi.h"
#include "openmc/constants.h"
#include "openmc/error.h"
#include "openmc/hdf5_interface.h"
#include "openmc/message_passing.h"
#include "openmc/search.h"
#include "openmc/settings.h"
#include "openmc/tallies/filter.h"
#include "openmc/xml_interface.h"

#ifdef LIBMESH
#include "libmesh/mesh_tools.h"
#endif

namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

namespace model {

std::vector<std::unique_ptr<Mesh>> meshes;
std::unordered_map<int32_t, int32_t> mesh_map;

} // namespace model

//==============================================================================
// Helper functions
//==============================================================================

//! Update an intersection point if the given candidate is closer.
//
//! The first 6 arguments are coordinates for the starting point of a particle
//! and its intersection with a mesh surface.  If the distance between these
//! two points is shorter than the given `min_distance`, then the `r` argument
//! will be updated to match the intersection point, and `min_distance` will
//! also be updated.

inline bool check_intersection_point(double x1, double x0, double y1,
  double y0, double z1, double z0, Position& r, double& min_distance)
{
  double dist = std::pow(x1-x0, 2) + std::pow(y1-y0, 2) + std::pow(z1-z0, 2);
  if (dist < min_distance) {
    r.x = x1;
    r.y = y1;
    r.z = z1;
    min_distance = dist;
    return true;
  }
  return false;
}

//==============================================================================
// Mesh implementation
//==============================================================================

Mesh::Mesh(pugi::xml_node node) {
   // Copy mesh id
  if (check_for_node(node, "id")) {
    id_ = std::stoi(get_node_value(node, "id"));

    // Check to make sure 'id' hasn't been used
    if (model::mesh_map.find(id_) != model::mesh_map.end()) {
      fatal_error("Two or more meshes use the same unique ID: " +
        std::to_string(id_));
    }
  }
}

//==============================================================================
// RegularMesh implementation
//==============================================================================

RegularMesh::RegularMesh(pugi::xml_node node)
  : Mesh {node}
{
  // Determine number of dimensions for mesh
  if (check_for_node(node, "dimension")) {
    shape_ = get_node_xarray<int>(node, "dimension");
    int n = n_dimension_ = shape_.size();
    if (n != 1 && n != 2 && n != 3) {
      fatal_error("Mesh must be one, two, or three dimensions.");
    }

    // Check that dimensions are all greater than zero
    if (xt::any(shape_ <= 0)) {
      fatal_error("All entries on the <dimension> element for a tally "
        "mesh must be positive.");
    }
  }

  // Check for lower-left coordinates
  if (check_for_node(node, "lower_left")) {
    // Read mesh lower-left corner location
    lower_left_ = get_node_xarray<double>(node, "lower_left");
  } else {
    fatal_error("Must specify <lower_left> on a mesh.");
  }

  if (check_for_node(node, "width")) {
    // Make sure both upper-right or width were specified
    if (check_for_node(node, "upper_right")) {
      fatal_error("Cannot specify both <upper_right> and <width> on a mesh.");
    }

    width_ = get_node_xarray<double>(node, "width");

    // Check to ensure width has same dimensions
    auto n = width_.size();
    if (n != lower_left_.size()) {
      fatal_error("Number of entries on <width> must be the same as "
        "the number of entries on <lower_left>.");
    }

    // Check for negative widths
    if (xt::any(width_ < 0.0)) {
      fatal_error("Cannot have a negative <width> on a tally mesh.");
    }

    // Set width and upper right coordinate
    upper_right_ = xt::eval(lower_left_ + shape_ * width_);

  } else if (check_for_node(node, "upper_right")) {
    upper_right_ = get_node_xarray<double>(node, "upper_right");

    // Check to ensure width has same dimensions
    auto n = upper_right_.size();
    if (n != lower_left_.size()) {
      fatal_error("Number of entries on <upper_right> must be the "
        "same as the number of entries on <lower_left>.");
    }

    // Check that upper-right is above lower-left
    if (xt::any(upper_right_ < lower_left_)) {
      fatal_error("The <upper_right> coordinates must be greater than "
        "the <lower_left> coordinates on a tally mesh.");
    }

    // Set width
    if (shape_.size() > 0) {
      width_ = xt::eval((upper_right_ - lower_left_) / shape_);
    }
  } else {
    fatal_error("Must specify either <upper_right> and <width> on a mesh.");
  }

  // Make sure lower_left and dimension match
  if (shape_.size() > 0) {
    if (shape_.size() != lower_left_.size()) {
      fatal_error("Number of entries on <lower_left> must be the same "
        "as the number of entries on <dimension>.");
    }

    // Set volume fraction
    volume_frac_ = 1.0/xt::prod(shape_)();
  }
}

int RegularMesh::get_bin(Position r) const
{
  // Loop over the dimensions of the mesh
  for (int i = 0; i < n_dimension_; ++i) {
    // Check for cases where particle is outside of mesh
    if (r[i] < lower_left_[i]) {
      return -1;
    } else if (r[i] > upper_right_[i]) {
      return -1;
    }
  }

  // Determine indices
  std::vector<int> ijk(n_dimension_);
  bool in_mesh;
  get_indices(r, ijk.data(), &in_mesh);
  if (!in_mesh) return -1;

  // Convert indices to bin
  return get_bin_from_indices(ijk.data());
}

int RegularMesh::get_bin_from_indices(const int* ijk) const
{
  switch (n_dimension_) {
  case 1:
    return ijk[0] - 1;
  case 2:
    return (ijk[1] - 1)*shape_[0] + ijk[0] - 1;
  case 3:
    return ((ijk[2] - 1)*shape_[1] + (ijk[1] - 1))*shape_[0] + ijk[0] - 1;
  default:
    throw std::runtime_error{"Invalid number of mesh dimensions"};
  }
}

void RegularMesh::get_indices(Position r, int* ijk, bool* in_mesh) const
{
  // Find particle in mesh
  *in_mesh = true;
  for (int i = 0; i < n_dimension_; ++i) {
    ijk[i] = std::ceil((r[i] - lower_left_[i]) / width_[i]);

    // Check if indices are within bounds
    if (ijk[i] < 1 || ijk[i] > shape_[i]) *in_mesh = false;
  }
}

void RegularMesh::get_indices_from_bin(int bin, int* ijk) const
{
  if (n_dimension_ == 1) {
    ijk[0] = bin + 1;
  } else if (n_dimension_ == 2) {
    ijk[0] = bin % shape_[0] + 1;
    ijk[1] = bin / shape_[0] + 1;
  } else if (n_dimension_ == 3) {
    ijk[0] = bin % shape_[0] + 1;
    ijk[1] = (bin % (shape_[0] * shape_[1])) / shape_[0] + 1;
    ijk[2] = bin / (shape_[0] * shape_[1]) + 1;
  }
}

int RegularMesh::n_bins() const
{
  int n_bins = 1;
  for (auto dim : shape_) n_bins *= dim;
  return n_bins;
}

int RegularMesh::n_surface_bins() const
{
  return 4 * n_dimension_ * n_bins();
}

bool RegularMesh::intersects(Position& r0, Position r1, int* ijk) const
{
  switch(n_dimension_) {
  case 1:
    return intersects_1d(r0, r1, ijk);
  case 2:
    return intersects_2d(r0, r1, ijk);
  case 3:
    return intersects_3d(r0, r1, ijk);
  default:
    throw std::runtime_error{"Invalid number of mesh dimensions."};
  }
}

bool RegularMesh::intersects_1d(Position& r0, Position r1, int* ijk) const
{
  // Copy coordinates of starting point
  double x0 = r0.x;
  double y0 = r0.y;
  double z0 = r0.z;

  // Copy coordinates of ending point
  double x1 = r1.x;
  double y1 = r1.y;
  double z1 = r1.z;

  // Copy coordinates of mesh lower_left and upper_right
  double xm0 = lower_left_[0];
  double xm1 = upper_right_[0];

  double min_dist = INFTY;

  // Check if line intersects left surface -- calculate the intersection point
  // (y,z)
  if ((x0 < xm0 && x1 > xm0) || (x0 > xm0 && x1 < xm0)) {
    double yi = y0 + (xm0 - x0) * (y1 - y0) / (x1 - x0);
    double zi = z0 + (xm0 - x0) * (z1 - z0) / (x1 - x0);
    if (check_intersection_point(xm0, x0, yi, yi, zi, zi, r0, min_dist)) {
      ijk[0] = 1;
    }
  }

  // Check if line intersects right surface -- calculate the intersection point
  // (y,z)
  if ((x0 < xm1 && x1 > xm1) || (x0 > xm1 && x1 < xm1)) {
    double yi = y0 + (xm1 - x0) * (y1 - y0) / (x1 - x0);
    double zi = z0 + (xm1 - x0) * (z1 - z0) / (x1 - x0);
    if (check_intersection_point(xm1, x0, yi, yi, zi, zi, r0, min_dist)) {
      ijk[0] = shape_[0];
    }
  }

  return min_dist < INFTY;
}

bool RegularMesh::intersects_2d(Position& r0, Position r1, int* ijk) const
{
  // Copy coordinates of starting point
  double x0 = r0.x;
  double y0 = r0.y;
  double z0 = r0.z;

  // Copy coordinates of ending point
  double x1 = r1.x;
  double y1 = r1.y;
  double z1 = r1.z;

  // Copy coordinates of mesh lower_left
  double xm0 = lower_left_[0];
  double ym0 = lower_left_[1];

  // Copy coordinates of mesh upper_right
  double xm1 = upper_right_[0];
  double ym1 = upper_right_[1];

  double min_dist = INFTY;

  // Check if line intersects left surface -- calculate the intersection point
  // (y,z)
  if ((x0 < xm0 && x1 > xm0) || (x0 > xm0 && x1 < xm0)) {
    double yi = y0 + (xm0 - x0) * (y1 - y0) / (x1 - x0);
    double zi = z0 + (xm0 - x0) * (z1 - z0) / (x1 - x0);
    if (yi >= ym0 && yi < ym1) {
      if (check_intersection_point(xm0, x0, yi, y0, zi, zi, r0, min_dist)) {
        ijk[0] = 1;
        ijk[1] = std::ceil((yi - lower_left_[1]) / width_[1]);
      }
    }
  }

  // Check if line intersects back surface -- calculate the intersection point
  // (x,z)
  if ((y0 < ym0 && y1 > ym0) || (y0 > ym0 && y1 < ym0)) {
    double xi = x0 + (ym0 - y0) * (x1 - x0) / (y1 - y0);
    double zi = z0 + (ym0 - y0) * (z1 - z0) / (y1 - y0);
    if (xi >= xm0 && xi < xm1) {
      if (check_intersection_point(xi, x0, ym0, y0, zi, zi, r0, min_dist)) {
        ijk[0] = std::ceil((xi - lower_left_[0]) / width_[0]);
        ijk[1] = 1;
      }
    }
  }

  // Check if line intersects right surface -- calculate the intersection point
  // (y,z)
  if ((x0 < xm1 && x1 > xm1) || (x0 > xm1 && x1 < xm1)) {
    double yi = y0 + (xm1 - x0) * (y1 - y0) / (x1 - x0);
    double zi = z0 + (xm1 - x0) * (z1 - z0) / (x1 - x0);
    if (yi >= ym0 && yi < ym1) {
      if (check_intersection_point(xm1, x0, yi, y0, zi, zi, r0, min_dist)) {
        ijk[0] = shape_[0];
        ijk[1] = std::ceil((yi - lower_left_[1]) / width_[1]);
      }
    }
  }

  // Check if line intersects front surface -- calculate the intersection point
  // (x,z)
  if ((y0 < ym1 && y1 > ym1) || (y0 > ym1 && y1 < ym1)) {
    double xi = x0 + (ym1 - y0) * (x1 - x0) / (y1 - y0);
    double zi = z0 + (ym1 - y0) * (z1 - z0) / (y1 - y0);
    if (xi >= xm0 && xi < xm1) {
      if (check_intersection_point(xi, x0, ym1, y0, zi, zi, r0, min_dist)) {
        ijk[0] = std::ceil((xi - lower_left_[0]) / width_[0]);
        ijk[1] = shape_[1];
      }
    }
  }

  return min_dist < INFTY;
}

bool RegularMesh::intersects_3d(Position& r0, Position r1, int* ijk) const
{
  // Copy coordinates of starting point
  double x0 = r0.x;
  double y0 = r0.y;
  double z0 = r0.z;

  // Copy coordinates of ending point
  double x1 = r1.x;
  double y1 = r1.y;
  double z1 = r1.z;

  // Copy coordinates of mesh lower_left
  double xm0 = lower_left_[0];
  double ym0 = lower_left_[1];
  double zm0 = lower_left_[2];

  // Copy coordinates of mesh upper_right
  double xm1 = upper_right_[0];
  double ym1 = upper_right_[1];
  double zm1 = upper_right_[2];

  double min_dist = INFTY;

  // Check if line intersects left surface -- calculate the intersection point
  // (y,z)
  if ((x0 < xm0 && x1 > xm0) || (x0 > xm0 && x1 < xm0)) {
    double yi = y0 + (xm0 - x0) * (y1 - y0) / (x1 - x0);
    double zi = z0 + (xm0 - x0) * (z1 - z0) / (x1 - x0);
    if (yi >= ym0 && yi < ym1 && zi >= zm0 && zi < zm1) {
      if (check_intersection_point(xm0, x0, yi, y0, zi, z0, r0, min_dist)) {
        ijk[0] = 1;
        ijk[1] = std::ceil((yi - lower_left_[1]) / width_[1]);
        ijk[2] = std::ceil((zi - lower_left_[2]) / width_[2]);
      }
    }
  }

  // Check if line intersects back surface -- calculate the intersection point
  // (x,z)
  if ((y0 < ym0 && y1 > ym0) || (y0 > ym0 && y1 < ym0)) {
    double xi = x0 + (ym0 - y0) * (x1 - x0) / (y1 - y0);
    double zi = z0 + (ym0 - y0) * (z1 - z0) / (y1 - y0);
    if (xi >= xm0 && xi < xm1 && zi >= zm0 && zi < zm1) {
      if (check_intersection_point(xi, x0, ym0, y0, zi, z0, r0, min_dist)) {
        ijk[0] = std::ceil((xi - lower_left_[0]) / width_[0]);
        ijk[1] = 1;
        ijk[2] = std::ceil((zi - lower_left_[2]) / width_[2]);
      }
    }
  }

  // Check if line intersects bottom surface -- calculate the intersection
  // point (x,y)
  if ((z0 < zm0 && z1 > zm0) || (z0 > zm0 && z1 < zm0)) {
    double xi = x0 + (zm0 - z0) * (x1 - x0) / (z1 - z0);
    double yi = y0 + (zm0 - z0) * (y1 - y0) / (z1 - z0);
    if (xi >= xm0 && xi < xm1 && yi >= ym0 && yi < ym1) {
      if (check_intersection_point(xi, x0, yi, y0, zm0, z0, r0, min_dist)) {
        ijk[0] = std::ceil((xi - lower_left_[0]) / width_[0]);
        ijk[1] = std::ceil((yi - lower_left_[1]) / width_[1]);
        ijk[2] = 1;
      }
    }
  }

  // Check if line intersects right surface -- calculate the intersection point
  // (y,z)
  if ((x0 < xm1 && x1 > xm1) || (x0 > xm1 && x1 < xm1)) {
    double yi = y0 + (xm1 - x0) * (y1 - y0) / (x1 - x0);
    double zi = z0 + (xm1 - x0) * (z1 - z0) / (x1 - x0);
    if (yi >= ym0 && yi < ym1 && zi >= zm0 && zi < zm1) {
      if (check_intersection_point(xm1, x0, yi, y0, zi, z0, r0, min_dist)) {
        ijk[0] = shape_[0];
        ijk[1] = std::ceil((yi - lower_left_[1]) / width_[1]);
        ijk[2] = std::ceil((zi - lower_left_[2]) / width_[2]);
      }
    }
  }

  // Check if line intersects front surface -- calculate the intersection point
  // (x,z)
  if ((y0 < ym1 && y1 > ym1) || (y0 > ym1 && y1 < ym1)) {
    double xi = x0 + (ym1 - y0) * (x1 - x0) / (y1 - y0);
    double zi = z0 + (ym1 - y0) * (z1 - z0) / (y1 - y0);
    if (xi >= xm0 && xi < xm1 && zi >= zm0 && zi < zm1) {
      if (check_intersection_point(xi, x0, ym1, y0, zi, z0, r0, min_dist)) {
        ijk[0] = std::ceil((xi - lower_left_[0]) / width_[0]);
        ijk[1] = shape_[1];
        ijk[2] = std::ceil((zi - lower_left_[2]) / width_[2]);
      }
    }
  }

  // Check if line intersects top surface -- calculate the intersection point
  // (x,y)
  if ((z0 < zm1 && z1 > zm1) || (z0 > zm1 && z1 < zm1)) {
    double xi = x0 + (zm1 - z0) * (x1 - x0) / (z1 - z0);
    double yi = y0 + (zm1 - z0) * (y1 - y0) / (z1 - z0);
    if (xi >= xm0 && xi < xm1 && yi >= ym0 && yi < ym1) {
      if (check_intersection_point(xi, x0, yi, y0, zm1, z0, r0, min_dist)) {
        ijk[0] = std::ceil((xi - lower_left_[0]) / width_[0]);
        ijk[1] = std::ceil((yi - lower_left_[1]) / width_[1]);
        ijk[2] = shape_[2];
      }
    }
  }

  return min_dist < INFTY;
}

void RegularMesh::bins_crossed(const Particle* p, std::vector<int>& bins,
                               std::vector<double>& lengths) const
{
  // ========================================================================
  // Determine where the track intersects the mesh and if it intersects at all.

  // Copy the starting and ending coordinates of the particle.
  Position last_r {p->r_last_};
  Position r {p->r()};
  Direction u {p->u()};

  // Compute the length of the entire track.
  double total_distance = (r - last_r).norm();

  // While determining if this track intersects the mesh, offset the starting
  // and ending coords by a bit.  This avoid finite-precision errors that can
  // occur when the mesh surfaces coincide with lattice or geometric surfaces.
  Position r0 = last_r + TINY_BIT*u;
  Position r1 = r - TINY_BIT*u;

  // Determine the mesh indices for the starting and ending coords.
  int n = n_dimension_;
  std::vector<int> ijk0(n), ijk1(n);
  bool start_in_mesh;
  get_indices(r0, ijk0.data(), &start_in_mesh);
  bool end_in_mesh;
  get_indices(r1, ijk1.data(), &end_in_mesh);

  // Reset coordinates and check for a mesh intersection if necessary.
  if (start_in_mesh) {
    // The initial coords lie in the mesh, use those coords for tallying.
    r0 = last_r;
  } else {
    // The initial coords do not lie in the mesh.  Check to see if the particle
    // eventually intersects the mesh and compute the relevant coords and
    // indices.
    if (!intersects(r0, r1, ijk0.data())) return;
  }
  r1 = r;

  // The TINY_BIT offsets above mean that the preceding logic cannot always find
  // the correct ijk0 and ijk1 indices. For tracks shorter than 2*TINY_BIT, just
  // assume the track lies in only one mesh bin. These tracks are very short so
  // any error caused by this assumption will be small.
  if (total_distance < 2*TINY_BIT) {
    for (int i = 0; i < n; ++i) ijk0[i] = ijk1[i];
  }

  // ========================================================================
  // Find which mesh cells are traversed and the length of each traversal.

  while (true) {
    if (ijk0 == ijk1) {
      // The track ends in this cell.  Use the particle end location rather
      // than the mesh surface and stop iterating.
      double distance = (r1 - r0).norm();
      bins.push_back(get_bin_from_indices(ijk0.data()));
      lengths.push_back(distance / total_distance);
      break;
    }

    // The track exits this cell.  Determine the distance to each mesh surface.
    std::vector<double> d(n);
    for (int k = 0; k < n; ++k) {
      if (std::fabs(u[k]) < FP_PRECISION) {
        d[k] = INFTY;
      } else if (u[k] > 0) {
        double xyz_cross = lower_left_[k] + ijk0[k] * width_[k];
        d[k] = (xyz_cross - r0[k]) / u[k];
      } else {
        double xyz_cross = lower_left_[k] + (ijk0[k] - 1) * width_[k];
        d[k] = (xyz_cross - r0[k]) / u[k];
      }
    }

    // Pick the closest mesh surface and append this traversal to the output.
    auto j = std::min_element(d.begin(), d.end()) - d.begin();
    double distance = d[j];
    bins.push_back(get_bin_from_indices(ijk0.data()));
    lengths.push_back(distance / total_distance);

    // Translate to the oncoming mesh surface.
    r0 += distance * u;

    // Increment the indices into the next mesh cell.
    if (u[j] > 0.0) {
      ++ijk0[j];
    } else {
      --ijk0[j];
    }

    // If the next indices are invalid, then the track has left the mesh and
    // we are done.
    bool in_mesh = true;
    for (int i = 0; i < n; ++i) {
      if (ijk0[i] < 1 || ijk0[i] > shape_[i]) {
        in_mesh = false;
        break;
      }
    }
    if (!in_mesh) break;
  }
}

void RegularMesh::surface_bins_crossed(const Particle* p,
                                       std::vector<int>& bins) const
{
  // ========================================================================
  // Determine if the track intersects the tally mesh.

  // Copy the starting and ending coordinates of the particle.
  Position r0 {p->r_last_current_};
  Position r1 {p->r()};
  Direction u {p->u()};

  // Determine indices for starting and ending location.
  int n = n_dimension_;
  std::vector<int> ijk0(n), ijk1(n);
  bool start_in_mesh;
  get_indices(r0, ijk0.data(), &start_in_mesh);
  bool end_in_mesh;
  get_indices(r1, ijk1.data(), &end_in_mesh);

  // Check if the track intersects any part of the mesh.
  if (!start_in_mesh) {
    Position r0_copy = r0;
    std::vector<int> ijk0_copy(ijk0);
    if (!intersects(r0_copy, r1, ijk0_copy.data())) return;
  }

  // ========================================================================
  // Find which mesh surfaces are crossed.

  // Calculate number of surface crossings
  int n_cross = 0;
  for (int i = 0; i < n; ++i) n_cross += std::abs(ijk1[i] - ijk0[i]);
  if (n_cross == 0) return;

  // Bounding coordinates
  Position xyz_cross;
  for (int i = 0; i < n; ++i) {
    if (u[i] > 0.0) {
      xyz_cross[i] = lower_left_[i] + ijk0[i] * width_[i];
    } else {
      xyz_cross[i] = lower_left_[i] + (ijk0[i] - 1) * width_[i];
    }
  }

  for (int j = 0; j < n_cross; ++j) {
    // Set the distances to infinity
    Position d {INFTY, INFTY, INFTY};

    // Determine closest bounding surface. We need to treat
    // special case where the cosine of the angle is zero since this would
    // result in a divide-by-zero.
    double distance = INFTY;
    for (int i = 0; i < n; ++i) {
      if (u[i] == 0) {
        d[i] = INFTY;
      } else {
        d[i] = (xyz_cross[i] - r0[i])/u[i];
      }
      distance = std::min(distance, d[i]);
    }

    // Loop over the dimensions
    for (int i = 0; i < n; ++i) {
      // Check whether distance is the shortest distance
      if (distance == d[i]) {

        // Check whether the current indices are within the mesh bounds
        bool in_mesh = true;
        for (int j = 0; j < n; ++j) {
          if (ijk0[j] < 1 || ijk0[j] > shape_[j]) {
            in_mesh = false;
            break;
          }
        }

        // Check whether particle is moving in positive i direction
        if (u[i] > 0) {

          // Outward current on i max surface
          if (in_mesh) {
            int i_surf = 4*i + 3;
            int i_mesh = get_bin_from_indices(ijk0.data());
            int i_bin = 4*n*i_mesh + i_surf - 1;

            bins.push_back(i_bin);
          }

          // Advance position
          ++ijk0[i];
          xyz_cross[i] += width_[i];
          in_mesh = true;
          for (int j = 0; j < n; ++j) {
            if (ijk0[j] < 1 || ijk0[j] > shape_[j]) {
              in_mesh = false;
              break;
            }
          }

          // If the particle crossed the surface, tally the inward current on
          // i min surface
          if (in_mesh) {
            int i_surf = 4*i + 2;
            int i_mesh = get_bin_from_indices(ijk0.data());
            int i_bin = 4*n*i_mesh + i_surf - 1;

            bins.push_back(i_bin);
          }

        } else {
          // The particle is moving in the negative i direction

          // Outward current on i min surface
          if (in_mesh) {
            int i_surf = 4*i + 1;
            int i_mesh = get_bin_from_indices(ijk0.data());
            int i_bin = 4*n*i_mesh + i_surf - 1;

            bins.push_back(i_bin);
          }

          // Advance position
          --ijk0[i];
          xyz_cross[i] -= width_[i];
          in_mesh = true;
          for (int j = 0; j < n; ++j) {
            if (ijk0[j] < 1 || ijk0[j] > shape_[j]) {
              in_mesh = false;
              break;
            }
          }

          // If the particle crossed the surface, tally the inward current on
          // i max surface
          if (in_mesh) {
            int i_surf = 4*i + 4;
            int i_mesh = get_bin_from_indices(ijk0.data());
            int i_bin = 4*n*i_mesh + i_surf - 1;

            bins.push_back(i_bin);
          }
        }
      }
    }

    // Calculate new coordinates
    r0 += distance * u;
  }
}

std::pair<std::vector<double>, std::vector<double>>
RegularMesh::plot(Position plot_ll, Position plot_ur) const
{
  // Figure out which axes lie in the plane of the plot.
  std::array<int, 2> axes {-1, -1};
  if (plot_ur.z == plot_ll.z) {
    axes[0] = 0;
    if (n_dimension_ > 1) axes[1] = 1;
  } else if (plot_ur.y == plot_ll.y) {
    axes[0] = 0;
    if (n_dimension_ > 2) axes[1] = 2;
  } else if (plot_ur.x == plot_ll.x) {
    if (n_dimension_ > 1) axes[0] = 1;
    if (n_dimension_ > 2) axes[1] = 2;
  } else {
    fatal_error("Can only plot mesh lines on an axis-aligned plot");
  }

  // Get the coordinates of the mesh lines along both of the axes.
  std::array<std::vector<double>, 2> axis_lines;
  for (int i_ax = 0; i_ax < 2; ++i_ax) {
    int axis = axes[i_ax];
    if (axis == -1) continue;
    auto& lines {axis_lines[i_ax]};

    double coord = lower_left_[axis];
    for (int i = 0; i < shape_[axis] + 1; ++i) {
      if (coord >= plot_ll[axis] && coord <= plot_ur[axis])
        lines.push_back(coord);
      coord += width_[axis];
    }
  }

  return {axis_lines[0], axis_lines[1]};
}

void RegularMesh::to_hdf5(hid_t group) const
{
  hid_t mesh_group = create_group(group, "mesh " + std::to_string(id_));

  write_dataset(mesh_group, "type", "regular");
  write_dataset(mesh_group, "dimension", shape_);
  write_dataset(mesh_group, "lower_left", lower_left_);
  write_dataset(mesh_group, "upper_right", upper_right_);
  write_dataset(mesh_group, "width", width_);

  close_group(mesh_group);
}

xt::xtensor<double, 1>
RegularMesh::count_sites(const std::vector<Particle::Bank>& bank,
                         //xt::xarray<double>
                         //Mesh::count_sites(const std::vector<Particle::Bank>& bank,
  bool* outside) const
{
  // Determine shape of array for counts
  std::size_t m = n_bins();
  std::vector<std::size_t> shape = {m};

  // Create array of zeros
  xt::xarray<double> cnt {shape, 0.0};
  bool outside_ = false;

  for (const auto& site : bank) {
    // determine scoring bin for entropy mesh
    int mesh_bin = get_bin(site.r);

    // if outside mesh, skip particle
    if (mesh_bin < 0) {
      outside_ = true;
      continue;
    }

    // Add to appropriate bin
    cnt(mesh_bin) += site.wgt;
  }

  // Create copy of count data. Since ownership will be acquired by xtensor,
  // std::allocator must be used to avoid Valgrind mismatched free() / delete
  // warnings.
  int total = cnt.size();
  double* cnt_reduced = std::allocator<double>{}.allocate(total);

#ifdef OPENMC_MPI
  // collect values from all processors
  MPI_Reduce(cnt.data(), cnt_reduced, total, MPI_DOUBLE, MPI_SUM, 0,
    mpi::intracomm);

  // Check if there were sites outside the mesh for any processor
  if (outside) {
    MPI_Reduce(&outside_, outside, 1, MPI_C_BOOL, MPI_LOR, 0, mpi::intracomm);
  }
#else
  std::copy(cnt.data(), cnt.data() + total, cnt_reduced);
  if (outside) *outside = outside_;
#endif

  // Adapt reduced values in array back into an xarray
  auto arr = xt::adapt(cnt_reduced, total, xt::acquire_ownership(), shape);
  xt::xarray<double> counts = arr;

  return counts;
}

//==============================================================================
// RectilinearMesh implementation
//==============================================================================

RectilinearMesh::RectilinearMesh(pugi::xml_node node)
  : Mesh {node}
{
  n_dimension_ = 3;

  grid_.resize(3);
  grid_[0] = get_node_array<double>(node, "x_grid");
  grid_[1] = get_node_array<double>(node, "y_grid");
  grid_[2] = get_node_array<double>(node, "z_grid");

  shape_ = {static_cast<int>(grid_[0].size()) - 1,
            static_cast<int>(grid_[1].size()) - 1,
            static_cast<int>(grid_[2].size()) - 1};

  for (const auto& g : grid_) {
    if (g.size() < 2) fatal_error("x-, y-, and z- grids for rectilinear meshes "
      "must each have at least 2 points");
    for (int i = 1; i < g.size(); ++i) {
      if (g[i] <= g[i-1]) fatal_error("Values in for x-, y-, and z- grids for "
        "rectilinear meshes must be sorted and unique.");
    }
  }

  lower_left_ = {grid_[0].front(), grid_[1].front(), grid_[2].front()};
  upper_right_ = {grid_[0].back(), grid_[1].back(), grid_[2].back()};
}

void RectilinearMesh::bins_crossed(const Particle* p, std::vector<int>& bins,
                                   std::vector<double>& lengths) const
{
  // ========================================================================
  // Determine where the track intersects the mesh and if it intersects at all.

  // Copy the starting and ending coordinates of the particle.
  Position last_r {p->r_last_};
  Position r {p->r()};
  Direction u {p->u()};

  // Compute the length of the entire track.
  double total_distance = (r - last_r).norm();

  // While determining if this track intersects the mesh, offset the starting
  // and ending coords by a bit.  This avoid finite-precision errors that can
  // occur when the mesh surfaces coincide with lattice or geometric surfaces.
  Position r0 = last_r + TINY_BIT*u;
  Position r1 = r - TINY_BIT*u;

  // Determine the mesh indices for the starting and ending coords.
  int ijk0[3], ijk1[3];
  bool start_in_mesh;
  get_indices(r0, ijk0, &start_in_mesh);
  bool end_in_mesh;
  get_indices(r1, ijk1, &end_in_mesh);

  // Reset coordinates and check for a mesh intersection if necessary.
  if (start_in_mesh) {
    // The initial coords lie in the mesh, use those coords for tallying.
    r0 = last_r;
  } else {
    // The initial coords do not lie in the mesh.  Check to see if the particle
    // eventually intersects the mesh and compute the relevant coords and
    // indices.
    if (!intersects(r0, r1, ijk0)) return;
  }
  r1 = r;

  // The TINY_BIT offsets above mean that the preceding logic cannot always find
  // the correct ijk0 and ijk1 indices. For tracks shorter than 2*TINY_BIT, just
  // assume the track lies in only one mesh bin. These tracks are very short so
  // any error caused by this assumption will be small.
  if (total_distance < 2*TINY_BIT) {
    for (int i = 0; i < 3; ++i) ijk0[i] = ijk1[i];
  }

  // ========================================================================
  // Find which mesh cells are traversed and the length of each traversal.

  while (true) {
    if (std::equal(ijk0, ijk0+3, ijk1)) {
      // The track ends in this cell.  Use the particle end location rather
      // than the mesh surface and stop iterating.
      double distance = (r1 - r0).norm();
      bins.push_back(get_bin_from_indices(ijk0));
      lengths.push_back(distance / total_distance);
      break;
    }

    // The track exits this cell.  Determine the distance to each mesh surface.
    double d[3];
    for (int k = 0; k < 3; ++k) {
      if (std::fabs(u[k]) < FP_PRECISION) {
        d[k] = INFTY;
      } else if (u[k] > 0) {
        double xyz_cross = grid_[k][ijk0[k]];
        d[k] = (xyz_cross - r0[k]) / u[k];
      } else {
        double xyz_cross = grid_[k][ijk0[k] - 1];
        d[k] = (xyz_cross - r0[k]) / u[k];
      }
    }

    // Pick the closest mesh surface and append this traversal to the output.
    auto j = std::min_element(d, d+3) - d;
    double distance = d[j];
    bins.push_back(get_bin_from_indices(ijk0));
    lengths.push_back(distance / total_distance);

    // Translate to the oncoming mesh surface.
    r0 += distance * u;

    // Increment the indices into the next mesh cell.
    if (u[j] > 0.0) {
      ++ijk0[j];
    } else {
      --ijk0[j];
    }

    // If the next indices are invalid, then the track has left the mesh and
    // we are done.
    bool in_mesh = true;
    for (int i = 0; i < 3; ++i) {
      if (ijk0[i] < 1 || ijk0[i] > shape_[i]) {
        in_mesh = false;
        break;
      }
    }
    if (!in_mesh) break;
  }
}

void RectilinearMesh::surface_bins_crossed(const Particle* p,
                                           std::vector<int>& bins) const
{
  // ========================================================================
  // Determine if the track intersects the tally mesh.

  // Copy the starting and ending coordinates of the particle.
  Position r0 {p->r_last_current_};
  Position r1 {p->r()};
  Direction u {p->u()};

  // Determine indices for starting and ending location.
  int ijk0[3], ijk1[3];
  bool start_in_mesh;
  get_indices(r0, ijk0, &start_in_mesh);
  bool end_in_mesh;
  get_indices(r1, ijk1, &end_in_mesh);

  // If the starting coordinates do not lie in the mesh, compute the coords and
  // mesh indices of the first intersection, and add the bin for this first
  // intersection.  Return if the particle does not intersect the mesh at all.
  if (!start_in_mesh) {
    // Compute the incoming intersection coordinates and indices.
    if (!intersects(r0, r1, ijk0)) return;

    // Determine which surface the particle entered.
    double min_dist = INFTY;
    int i_surf;
    for (int i = 0; i < 3; ++i) {
      if (u[i] > 0.0 && ijk0[i] == 1) {
        double d = std::abs(r0[i] - grid_[i][0]);
        if (d < min_dist) {
          min_dist = d;
          i_surf = 4*i + 2;
        }
      } else if (u[i] < 0.0 && ijk0[i] == shape_[i]) {
        double d = std::abs(r0[i] - grid_[i][shape_[i]]);
        if (d < min_dist) {
          min_dist = d;
          i_surf = 4*i + 4;
        }
      } // u[i] == 0 intentionally skipped
    }

    // Add the incoming current bin.
    int i_mesh = get_bin_from_indices(ijk0);
    int i_bin = 4*3*i_mesh + i_surf - 1;
    bins.push_back(i_bin);
  }

  // If the ending coordinates do not lie in the mesh, compute the coords and
  // mesh indices of the last intersection, and add the bin for this last
  // intersection.
  if (!end_in_mesh) {
    // Compute the outgoing intersection coordinates and indices.
    intersects(r1, r0, ijk1);

    // Determine which surface the particle exited.
    double min_dist = INFTY;
    int i_surf;
    for (int i = 0; i < 3; ++i) {
      if (u[i] > 0.0 && ijk1[i] == shape_[i]) {
        double d = std::abs(r1[i] - grid_[i][shape_[i]]);
        if (d < min_dist) {
          min_dist = d;
          i_surf = 4*i + 3;
        }
      } else if (u[i] < 0.0 && ijk1[i] == 1) {
        double d = std::abs(r1[i] - grid_[i][0]);
        if (d < min_dist) {
          min_dist = d;
          i_surf = 4*i + 1;
        }
      } // u[i] == 0 intentionally skipped
    }

    // Add the outgoing current bin.
    int i_mesh = get_bin_from_indices(ijk1);
    int i_bin = 4*3*i_mesh + i_surf - 1;
    bins.push_back(i_bin);
  }

  // ========================================================================
  // Find which mesh surfaces are crossed.

  // Calculate number of surface crossings
  int n_cross = 0;
  for (int i = 0; i < 3; ++i) n_cross += std::abs(ijk1[i] - ijk0[i]);
  if (n_cross == 0) return;

  // Bounding coordinates
  Position xyz_cross;
  for (int i = 0; i < 3; ++i) {
    if (u[i] > 0.0) {
      xyz_cross[i] = grid_[i][ijk0[i]];
    } else {
      xyz_cross[i] = grid_[i][ijk0[i] - 1];
    }
  }

  for (int j = 0; j < n_cross; ++j) {
    // Set the distances to infinity
    Position d {INFTY, INFTY, INFTY};

    // Determine closest bounding surface. We need to treat
    // special case where the cosine of the angle is zero since this would
    // result in a divide-by-zero.
    double distance = INFTY;
    for (int i = 0; i < 3; ++i) {
      if (u[i] == 0) {
        d[i] = INFTY;
      } else {
        d[i] = (xyz_cross[i] - r0[i])/u[i];
      }
      distance = std::min(distance, d[i]);
    }

    // Loop over the dimensions
    for (int i = 0; i < 3; ++i) {
      // Check whether distance is the shortest distance
      if (distance == d[i]) {

        // Check whether particle is moving in positive i direction
        if (u[i] > 0) {

          // Outward current on i max surface
          int i_surf = 4*i + 3;
          int i_mesh = get_bin_from_indices(ijk0);
          int i_bin = 4*3*i_mesh + i_surf - 1;
          bins.push_back(i_bin);

          // Advance position
          ++ijk0[i];
          xyz_cross[i] = grid_[i][ijk0[i]];

          // Inward current on i min surface
          i_surf = 4*i + 2;
          i_mesh = get_bin_from_indices(ijk0);
          i_bin = 4*3*i_mesh + i_surf - 1;
          bins.push_back(i_bin);

        } else {
          // The particle is moving in the negative i direction

          // Outward current on i min surface
          int i_surf = 4*i + 1;
          int i_mesh = get_bin_from_indices(ijk0);
          int i_bin = 4*3*i_mesh + i_surf - 1;
          bins.push_back(i_bin);

          // Advance position
          --ijk0[i];
          xyz_cross[i] = grid_[i][ijk0[i] - 1];

          // Inward current on i min surface
          i_surf = 4*i + 4;
          i_mesh = get_bin_from_indices(ijk0);
          i_bin = 4*3*i_mesh + i_surf - 1;
          bins.push_back(i_bin);
        }
      }
    }

    // Calculate new coordinates
    r0 += distance * u;
  }
}

int RectilinearMesh::get_bin(Position r) const
{
  // Determine indices
  int ijk[3];
  bool in_mesh;
  get_indices(r, ijk, &in_mesh);
  if (!in_mesh) return -1;

  // Convert indices to bin
  return get_bin_from_indices(ijk);
}

int RectilinearMesh::get_bin_from_indices(const int* ijk) const
{
  return ((ijk[2] - 1)*shape_[1] + (ijk[1] - 1))*shape_[0] + ijk[0] - 1;
}

void RectilinearMesh::get_indices(Position r, int* ijk, bool* in_mesh) const
{
  *in_mesh = true;

  for (int i = 0; i < 3; ++i) {
    if (r[i] < grid_[i].front() || r[i] > grid_[i].back()) {
      ijk[i] = -1;
      *in_mesh = false;
    } else {
      ijk[i] = lower_bound_index(grid_[i].begin(), grid_[i].end(), r[i]) + 1;
    }
  }
}

void RectilinearMesh::get_indices_from_bin(int bin, int* ijk) const
{
  ijk[0] = bin % shape_[0] + 1;
  ijk[1] = (bin % (shape_[0] * shape_[1])) / shape_[0] + 1;
  ijk[2] = bin / (shape_[0] * shape_[1]) + 1;
}

int RectilinearMesh::n_bins() const
{
  int n_bins = 1;
  for (auto dim : shape_) n_bins *= dim;
  return n_bins;
}

int RectilinearMesh::n_surface_bins() const
{
  return 4 * n_dimension_ * n_bins();
}

std::pair<std::vector<double>, std::vector<double>>
RectilinearMesh::plot(Position plot_ll, Position plot_ur) const
{
  // Figure out which axes lie in the plane of the plot.
  std::array<int, 2> axes {-1, -1};
  if (plot_ur.z == plot_ll.z) {
    axes = {0, 1};
  } else if (plot_ur.y == plot_ll.y) {
    axes = {0, 2};
  } else if (plot_ur.x == plot_ll.x) {
    axes = {1, 2};
  } else {
    fatal_error("Can only plot mesh lines on an axis-aligned plot");
  }

  // Get the coordinates of the mesh lines along both of the axes.
  std::array<std::vector<double>, 2> axis_lines;
  for (int i_ax = 0; i_ax < 2; ++i_ax) {
    int axis = axes[i_ax];
    std::vector<double>& lines {axis_lines[i_ax]};

    for (auto coord : grid_[axis]) {
      if (coord >= plot_ll[axis] && coord <= plot_ur[axis])
        lines.push_back(coord);
    }
  }

  return {axis_lines[0], axis_lines[1]};
}

void RectilinearMesh::to_hdf5(hid_t group) const
{
  hid_t mesh_group = create_group(group, "mesh " + std::to_string(id_));

  write_dataset(mesh_group, "type", "rectilinear");
  write_dataset(mesh_group, "x_grid", grid_[0]);
  write_dataset(mesh_group, "y_grid", grid_[1]);
  write_dataset(mesh_group, "z_grid", grid_[2]);

  close_group(mesh_group);
}

bool RectilinearMesh::intersects(Position& r0, Position r1, int* ijk) const
{
  // Copy coordinates of starting point
  double x0 = r0.x;
  double y0 = r0.y;
  double z0 = r0.z;

  // Copy coordinates of ending point
  double x1 = r1.x;
  double y1 = r1.y;
  double z1 = r1.z;

  // Copy coordinates of mesh lower_left
  double xm0 = grid_[0].front();
  double ym0 = grid_[1].front();
  double zm0 = grid_[2].front();

  // Copy coordinates of mesh upper_right
  double xm1 = grid_[0].back();
  double ym1 = grid_[1].back();
  double zm1 = grid_[2].back();

  double min_dist = INFTY;

  // Check if line intersects left surface -- calculate the intersection point
  // (y,z)
  if ((x0 < xm0 && x1 > xm0) || (x0 > xm0 && x1 < xm0)) {
    double yi = y0 + (xm0 - x0) * (y1 - y0) / (x1 - x0);
    double zi = z0 + (xm0 - x0) * (z1 - z0) / (x1 - x0);
    if (yi >= ym0 && yi < ym1 && zi >= zm0 && zi < zm1) {
      if (check_intersection_point(xm0, x0, yi, y0, zi, z0, r0, min_dist)) {
        ijk[0] = 1;
        ijk[1] = lower_bound_index(grid_[1].begin(), grid_[1].end(), yi) + 1;
        ijk[2] = lower_bound_index(grid_[2].begin(), grid_[2].end(), zi) + 1;
      }
    }
  }

  // Check if line intersects back surface -- calculate the intersection point
  // (x,z)
  if ((y0 < ym0 && y1 > ym0) || (y0 > ym0 && y1 < ym0)) {
    double xi = x0 + (ym0 - y0) * (x1 - x0) / (y1 - y0);
    double zi = z0 + (ym0 - y0) * (z1 - z0) / (y1 - y0);
    if (xi >= xm0 && xi < xm1 && zi >= zm0 && zi < zm1) {
      if (check_intersection_point(xi, x0, ym0, y0, zi, z0, r0, min_dist)) {
        ijk[0] = lower_bound_index(grid_[0].begin(), grid_[0].end(), xi) + 1;
        ijk[1] = 1;
        ijk[2] = lower_bound_index(grid_[2].begin(), grid_[2].end(), zi) + 1;
      }
    }
  }

  // Check if line intersects bottom surface -- calculate the intersection
  // point (x,y)
  if ((z0 < zm0 && z1 > zm0) || (z0 > zm0 && z1 < zm0)) {
    double xi = x0 + (zm0 - z0) * (x1 - x0) / (z1 - z0);
    double yi = y0 + (zm0 - z0) * (y1 - y0) / (z1 - z0);
    if (xi >= xm0 && xi < xm1 && yi >= ym0 && yi < ym1) {
      if (check_intersection_point(xi, x0, yi, y0, zm0, z0, r0, min_dist)) {
        ijk[0] = lower_bound_index(grid_[0].begin(), grid_[0].end(), xi) + 1;
        ijk[1] = lower_bound_index(grid_[1].begin(), grid_[1].end(), yi) + 1;
        ijk[2] = 1;
      }
    }
  }

  // Check if line intersects right surface -- calculate the intersection point
  // (y,z)
  if ((x0 < xm1 && x1 > xm1) || (x0 > xm1 && x1 < xm1)) {
    double yi = y0 + (xm1 - x0) * (y1 - y0) / (x1 - x0);
    double zi = z0 + (xm1 - x0) * (z1 - z0) / (x1 - x0);
    if (yi >= ym0 && yi < ym1 && zi >= zm0 && zi < zm1) {
      if (check_intersection_point(xm1, x0, yi, y0, zi, z0, r0, min_dist)) {
        ijk[0] = shape_[0];
        ijk[1] = lower_bound_index(grid_[1].begin(), grid_[1].end(), yi) + 1;
        ijk[2] = lower_bound_index(grid_[2].begin(), grid_[2].end(), zi) + 1;
      }
    }
  }

  // Check if line intersects front surface -- calculate the intersection point
  // (x,z)
  if ((y0 < ym1 && y1 > ym1) || (y0 > ym1 && y1 < ym1)) {
    double xi = x0 + (ym1 - y0) * (x1 - x0) / (y1 - y0);
    double zi = z0 + (ym1 - y0) * (z1 - z0) / (y1 - y0);
    if (xi >= xm0 && xi < xm1 && zi >= zm0 && zi < zm1) {
      if (check_intersection_point(xi, x0, ym1, y0, zi, z0, r0, min_dist)) {
        ijk[0] = lower_bound_index(grid_[0].begin(), grid_[0].end(), xi) + 1;
        ijk[1] = shape_[1];
        ijk[2] = lower_bound_index(grid_[2].begin(), grid_[2].end(), zi) + 1;
      }
    }
  }

  // Check if line intersects top surface -- calculate the intersection point
  // (x,y)
  if ((z0 < zm1 && z1 > zm1) || (z0 > zm1 && z1 < zm1)) {
    double xi = x0 + (zm1 - z0) * (x1 - x0) / (z1 - z0);
    double yi = y0 + (zm1 - z0) * (y1 - y0) / (z1 - z0);
    if (xi >= xm0 && xi < xm1 && yi >= ym0 && yi < ym1) {
      if (check_intersection_point(xi, x0, yi, y0, zm1, z0, r0, min_dist)) {
        ijk[0] = lower_bound_index(grid_[0].begin(), grid_[0].end(), xi) + 1;
        ijk[1] = lower_bound_index(grid_[1].begin(), grid_[1].end(), yi) + 1;
        ijk[2] = shape_[2];
      }
    }
  }

  return min_dist < INFTY;
}

//==============================================================================
// Helper functions for the C API
//==============================================================================

int
check_mesh(int32_t index)
{
  if (index < 0 || index >= model::meshes.size()) {
    set_errmsg("Index in meshes array is out of bounds.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }
  return 0;
}

int
check_regular_mesh(int32_t index, RegularMesh** mesh)
{
  if (int err = check_mesh(index)) return err;
  *mesh = dynamic_cast<RegularMesh*>(model::meshes[index].get());
  if (!*mesh) {
    set_errmsg("This function is only valid for regular meshes.");
    return OPENMC_E_INVALID_TYPE;
  }
  return 0;
}

//==============================================================================
// C API functions
//==============================================================================

RegularMesh* get_regular_mesh(int32_t index) {
  return dynamic_cast<RegularMesh*>(model::meshes[index].get());
}

//! Extend the meshes array by n elements
extern "C" int
openmc_extend_meshes(int32_t n, int32_t* index_start, int32_t* index_end)
{
  if (index_start) *index_start = model::meshes.size();
  for (int i = 0; i < n; ++i) {
    model::meshes.push_back(std::move(std::make_unique<RegularMesh>()));
  }
  if (index_end) *index_end = model::meshes.size() - 1;

  return 0;
}

//! Return the index in the meshes array of a mesh with a given ID
extern "C" int
openmc_get_mesh_index(int32_t id, int32_t* index)
{
  auto pair = model::mesh_map.find(id);
  if (pair == model::mesh_map.end()) {
    set_errmsg("No mesh exists with ID=" + std::to_string(id) + ".");
    return OPENMC_E_INVALID_ID;
  }
  *index = pair->second;
  return 0;
}

// Return the ID of a mesh
extern "C" int
openmc_mesh_get_id(int32_t index, int32_t* id)
{
  if (int err = check_mesh(index)) return err;
  *id = model::meshes[index]->id_;
  return 0;
}

//! Set the ID of a mesh
extern "C" int
openmc_mesh_set_id(int32_t index, int32_t id)
{
  if (int err = check_mesh(index)) return err;
  model::meshes[index]->id_ = id;
  model::mesh_map[id] = index;
  return 0;
}

//! Get the dimension of a mesh
extern "C" int
openmc_mesh_get_dimension(int32_t index, int** dims, int* n)
{
  if (index < 0 || index >= model::meshes.size()) {
    set_errmsg("Index in meshes array is out of bounds.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }
  RegularMesh* mesh;
  if (int err = check_regular_mesh(index, &mesh)) return err;
  *dims = mesh->shape_.data();
  *n = mesh->n_dimension_;
  return 0;
}

//! Set the dimension of a mesh
extern "C" int
openmc_mesh_set_dimension(int32_t index, int n, const int* dims)
{
  if (index < 0 || index >= model::meshes.size()) {
    set_errmsg("Index in meshes array is out of bounds.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }

  RegularMesh* mesh;
  if (int err = check_regular_mesh(index, &mesh)) return err;

  // Copy dimension
  std::vector<std::size_t> shape = {static_cast<std::size_t>(n)};
  mesh->shape_ = xt::adapt(dims, n, xt::no_ownership(), shape);
  mesh->n_dimension_ = mesh->shape_.size();
  return 0;
}

//! Get the mesh parameters
extern "C" int
openmc_mesh_get_params(int32_t index, double** ll, double** ur, double** width, int* n)
{
  RegularMesh* m;
  if (int err = check_regular_mesh(index, &m)) return err;

  if (m->lower_left_.dimension() == 0) {
    set_errmsg("Mesh parameters have not been set.");
    return OPENMC_E_ALLOCATE;
  }

  *ll = m->lower_left_.data();
  *ur = m->upper_right_.data();
  *width = m->width_.data();
  *n = m->n_dimension_;
  return 0;
}

//! Set the mesh parameters
extern "C" int
openmc_mesh_set_params(int32_t index, int n, const double* ll, const double* ur,
                       const double* width)
{
  RegularMesh* m;
  if (int err = check_regular_mesh(index, &m)) return err;

  std::vector<std::size_t> shape = {static_cast<std::size_t>(n)};
  if (ll && ur) {
    m->lower_left_ = xt::adapt(ll, n, xt::no_ownership(), shape);
    m->upper_right_ = xt::adapt(ur, n, xt::no_ownership(), shape);
    m->width_ = (m->upper_right_ - m->lower_left_) / m->shape_;
  } else if (ll && width) {
    m->lower_left_ = xt::adapt(ll, n, xt::no_ownership(), shape);
    m->width_ = xt::adapt(width, n, xt::no_ownership(), shape);
    m->upper_right_ = m->lower_left_ + m->shape_ * m->width_;
  } else if (ur && width) {
    m->upper_right_ = xt::adapt(ur, n, xt::no_ownership(), shape);
    m->width_ = xt::adapt(width, n, xt::no_ownership(), shape);
    m->lower_left_ = m->upper_right_ - m->shape_ * m->width_;
  } else {
    set_errmsg("At least two parameters must be specified.");
    return OPENMC_E_INVALID_ARGUMENT;
  }

  return 0;
}

#ifdef DAGMC

UnstructuredMesh::UnstructuredMesh(pugi::xml_node node) : Mesh(node) {
  // check the mesh type
  if (check_for_node(node, "type")) {
    auto temp = get_node_value(node, "type", true, true);
    if (temp != "unstructured") {
      fatal_error("Invalid mesh type: " + temp);
    }
  }

  // get the filename of the unstructured mesh to load
  if (check_for_node(node, "mesh_file")) {
    filename_ = get_node_value(node, "mesh_file");
  }
  else {
    fatal_error("No filename supplied for unstructured mesh with ID: " +
                std::to_string(id_));
  }

  // create MOAB instance
  mbi_ = std::shared_ptr<moab::Interface>(new moab::Core());
  // create meshset to load mesh into
  moab::ErrorCode rval = mbi_->create_meshset(moab::MESHSET_SET, meshset_);
  if (rval != moab::MB_SUCCESS) {
    fatal_error("Failed to create fileset for umesh: " + filename_);
  }
  // load unstructured mesh file
  rval = mbi_->load_file(filename_.c_str(), &meshset_);
  if (rval != moab::MB_SUCCESS) {
    fatal_error("Failed to load the unstructured mesh file: " + filename_);
  }

  // always 3 for unstructured meshes
  n_dimension_ = 3;

  moab::Range all_tets;
  rval = mbi_->get_entities_by_dimension(meshset_, n_dimension_, all_tets);
  if (rval != moab::MB_SUCCESS) {
    fatal_error("Failed to get all tetrahedral elements");
  }

  ehs_.clear();
  ehs_ = all_tets;

  if (!all_tets.all_of_type(moab::MBTET)) {
    warning("Non-tetrahedral elements found in unstructured mesh: " + filename_);
  }

  compute_barycentric_data(ehs_);
  build_kdtree(ehs_);
}

void
UnstructuredMesh::build_kdtree(const moab::Range& all_tets) {

  moab::Range all_tris;
  int triangle_dim = 2;
  moab::ErrorCode rval = mbi_->get_adjacencies(all_tets,
                                               triangle_dim,
                                               true,
                                               all_tris,
                                               moab::Interface::UNION);
  if (rval != moab::MB_SUCCESS) {
    fatal_error("Failed to get adjacent triangles for tets");
  }

  if (!all_tris.all_of_type(moab::MBTRI)) {
    warning("Non-triangle elements found in tet adjacencies in unstructured mesh: " + filename_);
  }

  // combine into one range
  moab::Range all_tets_and_tris;
  all_tets_and_tris.merge(all_tets);
  all_tets_and_tris.merge(all_tris);

  // create and build KD-tree
  kdtree_ = std::unique_ptr<moab::AdaptiveKDTree>(new moab::AdaptiveKDTree(mbi_.get()));

  // build the tree
  rval = kdtree_->build_tree(all_tets_and_tris, &kdtree_root_);
  if (rval != moab::MB_SUCCESS) {
    fatal_error("Failed to construct KDTree for the unstructured mesh " + filename_);
  }
}

void
UnstructuredMesh::intersect_track(const moab::CartVect& start,
                                  const moab::CartVect& dir,
                                  double track_len,
                                  UnstructuredMeshHits& hits) const {
  moab::ErrorCode rval;
  std::vector<moab::EntityHandle> tris;
  std::vector<double> intersection_dists;

  rval = kdtree_->ray_intersect_triangles(kdtree_root_,
                                          1E-06,
                                          dir.array(),
                                          start.array(),
                                          tris,
                                          intersection_dists,
                                          0,
                                          track_len);
  if (rval != moab::MB_SUCCESS) {
    fatal_error("Failed to compute intersections on umesh: " + filename_);
  }

  // sort the tris and intersections by distance
  hits.clear();
  for (int i = 0; i < tris.size(); i++) {
    hits.push_back(std::pair<double, moab::EntityHandle>(intersection_dists[i], tris[i]));
  }

  // sorts by first component of std::pair by default
  std::sort(hits.begin(), hits.end());

  for (const auto& hit  : hits) {

  }

  for (const auto& hit : hits) {
    moab::CartVect hit_loc = start + dir * hit.first;
  }

}

void
UnstructuredMesh::bins_crossed(const Particle* p,
                               std::vector<int>& bins,
                               std::vector<double>& lengths) const {
  moab::ErrorCode rval;
  Position last_r{p->r_last_};
  Position r{p->r()};
  Position u{p->u()};
  u /= u.norm();
  moab::CartVect r0(last_r.x, last_r.y, last_r.z);
  moab::CartVect r1(r.x, r.y, r.z);
  moab::CartVect dir(u.x, u.y, u.z);

  double track_len = (r1 - r0).length();

  r0 += TINY_BIT*dir;
  r1 -= TINY_BIT*dir;

  UnstructuredMeshHits hits;
  intersect_track(r0, dir, track_len, hits);

  bins.clear();
  lengths.clear();

  // track could be entirely contained by a tet
  if (hits.size() == 0) {
    auto last_r_tet = get_tet(last_r + u * track_len * 0.5);
    if (last_r_tet) {
      bins.push_back(get_bin_from_ent_handle(last_r_tet));
      lengths.push_back(1.0);
    }
    return;
  }

  double last_dist = 0.0;
  for (auto hit = hits.begin(); hit != hits.end(); hit++) {

    // mid point for this segment
    double segment_length = hit->first - last_dist;
    Position segment_midpoint = last_r + u * last_dist + u * segment_length / 2.0;

    last_dist = hit->first;

    // try to get a tet at the midpoint of this segment
    auto tet = get_tet(segment_midpoint);
    if (tet) {
      bins.push_back(get_bin_from_ent_handle(tet));
      lengths.push_back(segment_length / track_len);
    }
  }

  // tally remaining portion of track after last hit if
  // the last segment of the track is in the mesh
  if (hits.back().first < track_len) {
    auto pos = (last_r + u * hits.back().first) + u * ((track_len - hits.back().first) / 2.0);
    auto tet = get_tet(pos);
    if (tet) {
      bins.push_back(get_bin_from_ent_handle(tet));
      lengths.push_back((track_len - hits.back().first) / track_len);
    }
  }

  return;
}

moab::EntityHandle
UnstructuredMesh::get_tet(const Position& r) const {
  moab::CartVect pos(r.x, r.y, r.z);
  moab::AdaptiveKDTreeIter kdtree_iter;
  moab::ErrorCode rval = kdtree_->point_search(pos.array(), kdtree_iter);
  if (rval != moab::MB_SUCCESS) { return 0; }

  moab::EntityHandle leaf = kdtree_iter.handle();
  moab::Range tets;
  rval = mbi_->get_entities_by_dimension(leaf, 3, tets, false);
  if (rval != moab::MB_SUCCESS) {
    warning("MOAB error finding tets.");
  }

  for (const auto& tet : tets) {
      if (point_in_tet(pos, tet)) {
        return tet;
    }
  }
  return 0;
}

double UnstructuredMesh::tet_volume(moab::EntityHandle tet) const {
 std::vector<moab::EntityHandle> conn;
 moab::ErrorCode rval = mbi_->get_connectivity(&tet, 1, conn);
 if (rval != moab::MB_SUCCESS) {
   fatal_error("Failed to get tet connectivity");
 }

 moab::CartVect p[4];
 rval = mbi_->get_coords(&(conn[0]), (int)conn.size(), p[0].array());
 if (rval != moab::MB_SUCCESS) {
   fatal_error("Failed to get tet coords");
 }

 return 1.0 / 6.0 * (((p[1] - p[0]) * (p[2] - p[0])) % (p[3] - p[0]));
}

//! Determine which surface bins were crossed by a particle
//
//! \param[in] p Particle to check
//! \param[out] bins Surface bins that were crossed
void UnstructuredMesh::surface_bins_crossed(const Particle* p, std::vector<int>& bins) const {
  return;
}

std::string UnstructuredMesh::get_label_for_bin(int bin) const {
  std::stringstream out;
  out << "MOAB EntityHandle: " << get_ent_handle_from_bin(bin);
  return out.str();
}

int
UnstructuredMesh::get_bin(Position r) const {
  moab::EntityHandle tet = get_tet(r);
  if (tet == 0) {
    return -1;
  } else {
    return get_bin_from_ent_handle(tet);
  }
}

void
UnstructuredMesh::compute_barycentric_data(const moab::Range& all_tets) {
  moab::ErrorCode rval;

  baryc_data_.clear();
  baryc_data_.resize(all_tets.size());

  for (auto& tet : all_tets) {
    std::vector<moab::EntityHandle> verts;
    rval = mbi_->get_connectivity(&tet, 1, verts);
    if (rval != moab::MB_SUCCESS) {
      fatal_error("Failed to get connectivity of tet on umesh: " + filename_);
    }

    moab::CartVect p[4];
    rval = mbi_->get_coords(&(verts[0]), (int)verts.size(), p[0].array());
    if (rval != moab::MB_SUCCESS) {
      fatal_error("Failed to get coordinates of a tet in umesh: " + filename_);
    }

    moab::Matrix3 a(p[1] - p[0], p[2] - p[0], p[3] - p[0], true);
    a = a.transpose().inverse();
    baryc_data_.at(get_bin_from_ent_handle(tet)) = a;
  }

}

// TODO: write this function
void
UnstructuredMesh::to_hdf5(hid_t group) const
{
    hid_t mesh_group = create_group(group, "mesh " + std::to_string(id_));

    write_dataset(mesh_group, "type", "unstructured");
    write_dataset(mesh_group, "filename", filename_);
    write_dataset(mesh_group, "library", "moab");
    // write volume of each tet
    std::vector<double> tet_vols;
    for (const auto& eh : ehs_) {
      tet_vols.emplace_back(tet_volume(eh));
    }
    write_dataset(mesh_group, "volumes", tet_vols);

    close_group(mesh_group);
}

bool
UnstructuredMesh::point_in_tet(const moab::CartVect& r, moab::EntityHandle tet) const {

  moab::ErrorCode rval;

  // get tet vertices
  std::vector<moab::EntityHandle> verts;
  rval = mbi_->get_connectivity(&tet, 1, verts);
  if (rval != moab::MB_SUCCESS) {
    warning("Failed to get vertices of tet in umesh: " + filename_);
    return false;
  }

  moab::EntityHandle v_zero = verts[0];
  moab::CartVect p_zero;
  rval = mbi_->get_coords(&v_zero, 1, p_zero.array());
  if (rval != moab::MB_SUCCESS) {
    warning("Failed to get coordinates of a vertex in umesh: " + filename_);
    return false;
  }

  // look up barycentric data
  int idx = get_bin_from_ent_handle(tet);
  const moab::Matrix3& a_inv = baryc_data_[idx];

  moab::CartVect bary_coords = a_inv * (r - p_zero);

  bool in_tet = (bary_coords[0] >= 0 && bary_coords[1] >= 0 && bary_coords[2] >= 0 &&
                 bary_coords[0] + bary_coords[1] + bary_coords[2] <= 1.);

  return in_tet;
}

int
UnstructuredMesh::get_bin_from_indices(const int* ijk) const {
  if (ijk[0] >= n_bins()) {
    std::stringstream s;
    s << "Invalid bin: " << ijk[0];
    fatal_error(s);
  }
  int bin = ehs_[ijk[0]] - ehs_[0];
  return bin;
}

void
UnstructuredMesh::get_indices(Position r, int* ijk, bool* in_mesh) const {
  int bin = get_bin(r);
  ijk[0]= bin;
  *in_mesh = bin != -1;
}

void UnstructuredMesh::get_indices_from_bin(int bin, int* ijk) const {
  ijk[0] = bin;
}

std::pair<std::vector<double>, std::vector<double>>
UnstructuredMesh::plot(Position plot_ll, Position plot_ur) const { return {}; }

int
UnstructuredMesh::get_bin_from_ent_handle(moab::EntityHandle eh) const {
  int bin = eh - ehs_[0];
  if (bin >= n_bins()) {
    std::stringstream s;
    s << "Invalid bin: " << bin;
    fatal_error(s);
  }
  return bin;
}

moab::EntityHandle
UnstructuredMesh::get_ent_handle_from_bin(int bin) const {
  if (bin >= n_bins()) {
    std::stringstream s;
    s << "Invalid bin index: " << bin;
    fatal_error(s);
  }
  return ehs_[0] + bin;
}

int UnstructuredMesh::n_bins() const {
  return ehs_.size();
}

int UnstructuredMesh::n_surface_bins() const {
  // collect all triangles in the set of tets for this mesh
  moab::Range tris;
  moab::ErrorCode rval;
  rval = mbi_->get_entities_by_type(0, moab::MBTRI, tris);
  if (rval != moab::MB_SUCCESS) {
    warning("Failed to get all triangles in the mesh instance");
    return -1;
  }
  return 2 * tris.size();
}

#endif

UnstructuredMeshBase::UnstructuredMeshBase(pugi::xml_node node) : Mesh(node) {
    // check the mesh type
  if (check_for_node(node, "type")) {
    auto temp = get_node_value(node, "type", true, true);
    if (temp != "unstructured") {
      fatal_error("Invalid mesh type: " + temp);
    }
  }

  // get the filename of the unstructured mesh to load
  if (check_for_node(node, "mesh_file")) {
    filename_ = get_node_value(node, "mesh_file");
  }
  else {
    fatal_error("No filename supplied for unstructured mesh with ID: " +
                std::to_string(id_));
  }
}

#ifdef LIBMESH
LibMesh::LibMesh(pugi::xml_node node) : UnstructuredMeshBase(node) {

  // always 3 for unstructured meshes
  n_dimension_ = 3;

  m_ = std::unique_ptr<libMesh::Mesh>(new libMesh::Mesh(settings::LMI->comm(), 3));
  m_->read(filename_);
  m_->prepare_for_use();

  eq_system_name_ = "mesh_" + std::to_string(id_) + "_system";

  equation_systems_ =
    std::unique_ptr<libMesh::EquationSystems>(new libMesh::EquationSystems(*m_));
  libMesh::ExplicitSystem& eq_sys =
    equation_systems_->add_system<libMesh::ExplicitSystem>(eq_system_name_);

  // setup the point locator
  int n_threads = omp_get_max_threads();
  point_locators_ = std::vector<std::unique_ptr<libMesh::PointLocatorBase>>(n_threads);
  for (int i = 0; i < n_threads; i++) {
    point_locators_[i] = m_->sub_point_locator();
    point_locators_[i]->set_find_element_tol(FP_COINCIDENT);
    point_locators_[i]->enable_out_of_mesh_mode();
  }

  // will need mesh neighbors to walk the mesh
  m_->find_neighbors();

  first_element_ = *m_->elements_begin();

  // bounding box for the mesh
  bbox_ = libMesh::MeshTools::create_bounding_box(*m_);
  bsphere_ = libMesh::MeshTools::bounding_sphere(*m_);

  // determine boundary elements and create bounding box
  for (int i = 0; i < m_->n_elem(); i++) {
    auto e = m_->elem_ptr(i);

    for (int k = 0; k < e->n_neighbors(); k++) {
      if (!e->neighbor_ptr(k)) {
        boundary_elements_.insert(e);
      }
    }
  }
}

void
LibMesh::get_indices(Position r, int* ijk, bool* in_mesh) const {
  int bin = get_bin(r);
  *in_mesh = bin != -1;
  ijk[0] = bin;
  return;
}

int LibMesh::n_bins() const {
  return m_->n_elem();
}

int LibMesh::n_surface_bins() const {
  return 0;
}

void
LibMesh::surface_bins_crossed(const Particle* p,
                               std::vector<int>& bins) const
{}

std::pair<std::vector<double>, std::vector<double>>
LibMesh::plot(Position plot_ll,
              Position plot_ur) const { return {}; }


int
LibMesh::get_bin_from_indices(const int* ijk) const {
  if (ijk[0] >= n_bins()) {
    std::stringstream s;
    s << "Invalid bin: " << ijk[0];
    fatal_error(s);
  }
  int bin = first_element_->id() + ijk[0];
  return bin;
}

void
LibMesh::get_indices_from_bin(int bin, int* ijk) const {
  ijk[0] = bin;
}

const libMesh::Elem*
LibMesh::get_element_from_bin(int bin) const {
  m_->elem_ptr(bin);
}

void
LibMesh::add_variable(const std::string& var_name) {
  // check if this is a new varaible
  if (!variable_map_.count(var_name)) {
    auto& eqn_sys = equation_systems_->get_system(eq_system_name_);
    auto var_num = eqn_sys.add_variable(var_name, libMesh::CONSTANT, libMesh::MONOMIAL);
    variable_map_[var_name] = var_num;
    equation_systems_->init();
  }
}

void
LibMesh::set_variable(const std::string& var_name, int bin, double value) {
  // look up the variable
  unsigned int var_num = variable_map_.at(var_name);
  auto& eqn_sys = equation_systems_->get_system(eq_system_name_);
  const libMesh::DofMap& dof_map = eqn_sys.get_dof_map();
  auto e = m_->elem_ptr(bin);
  std::vector<libMesh::dof_id_type> dof_indices;
  dof_map.dof_indices(e, dof_indices, var_num);
  Ensures(dof_indices.size() == 1);
  (*eqn_sys.solution).set(dof_indices[0], value);
}

void LibMesh::write(const std::string& filename) const {
  libMesh::ExodusII_IO exo(*m_);
  std::set<std::string> systems_out = {eq_system_name_};
  exo.write_discontinuous_exodusII(filename, *equation_systems_, &systems_out);
}


void
LibMesh::bins_crossed(const Particle* p,
                      std::vector<int>& bins,
                      std::vector<double>& lengths) const
{
  // get element containing previous position
  libMesh::Point start(p->r_last_.x, p->r_last_.y, p->r_last_.z);
  libMesh::Point end(p->r().x, p->r().y, p->r().z);
  libMesh::Point dir(p->u().x, p->u().y, p->u().z);
  dir /= dir.norm();

  double track_len = (end - start).norm();

  UnstructuredMeshHits hits;
  intersect_track(start, dir, track_len, hits);

  bins.clear();
  lengths.clear();

  for (const auto& hit : hits) {
    lengths.push_back(hit.first / track_len);
    bins.push_back(get_bin_from_element(hit.second));
  }
}

int
LibMesh::get_bin(Position r) const
{
  // look-up a tet using the point locator
  libMesh::Point p(r.x, r.y, r.z);

  // quick rejection check
  //  if (bsphere_.above_surface(p) || !bbox_.contains_point(p)) { return -1; }
  if (!bbox_.contains_point(p)) { return -1; }
  int thread = omp_get_thread_num();
  auto e = (*point_locators_[thread])(p);

  if (!e) {
    return -1;
  } else {
    return get_bin_from_element(e);
  }
}

bool LibMesh::intersects(Position& r0, Position r1, int* ijk) const {

  // first try to locate an element
  // for the start or end point
  int bin {-1};
  bin = get_bin(r0);
  if (bin != -1) {
    ijk[0] = bin;
    return true;
  }

  bin = get_bin(r1);
  if (bin != -1) {
    ijk[0] = bin;
    return true;
  }

  // check for an intersection with one of the boundary faces
  auto result = locate_boundary_element(r0, r1);
  // if we don't get a hit, the track won't intersect with the mesh
  if (result.second) {
    ijk[0] = get_bin_from_element(result.second);
    return true;
  }

  return false;
}

std::pair<double, const libMesh::Elem*>
LibMesh::locate_boundary_element(const Position& r0,
                                 const Position& r1) const {
  libMesh::Point a(r0.x, r0.y, r0.z);
  libMesh::Point b(r1.x, r1.y, r1.z);

  return locate_boundary_element(a, b);
}

std::pair<double, const libMesh::Elem*>
LibMesh::locate_boundary_element(const libMesh::Point& start,
                                 const libMesh::Point& end) const
{
  typedef std::pair<double, const libMesh::Elem*> RayHit;
  RayHit result = {INFTY, nullptr};
  if (start == end) { return result; }

  // attempt to locate an intersection with the mesh boundary
  libMesh::Point dir = (end - start).unit();
  double length = (end - start).norm();

  // locate potential elements
  std::set<const libMesh::Elem*> candidate_elements;
  for (auto elem : boundary_elements_) {
    // being conservative about search parameter by adding element hmax
    if (elem->close_to_point(start, length + elem->hmax())) {
      candidate_elements.insert(elem);
    }
  }

  // find nearest hit along our direction
  for (auto elem : candidate_elements) {
    for (int i = 0; i < elem->n_sides(); i++) {
      double temp_dist = 0;
      bool hit = plucker_test(elem->side_ptr(i), start, dir, temp_dist);
      if (hit && temp_dist > FP_COINCIDENT && temp_dist <= length) {
        // update if we find a closer intersection
        if (temp_dist < result.first) { result = {temp_dist, elem}; }
      }
    }
  }

  return result;
}

int
LibMesh::get_bin_from_element(const libMesh::Elem* elem) const {
  int bin =  elem->id() - first_element_->id();
  if (bin >= n_bins()) {
    std::stringstream s;
    s << "Invalid bin: " << bin;
    fatal_error(s);
  }
  return bin;
}

bool
LibMesh::inside_tet(const libMesh::Point& r,
                    const libMesh::Point& u,
                    std::unique_ptr<libMesh::Elem> e) const
{
  return inside_tet(r, u, e.get());
}


bool
LibMesh::inside_tet(const libMesh::Point& r,
                    const libMesh::Point& u,
                    const libMesh::Elem* e) const
{
  // fire rays at each triangle in the tet
  int n_hits = 0;
  for (int i = 0; i < e->n_sides(); i++) {
    double temp;
    if (plucker_test(e->side_ptr(i), r, u, temp)) { n_hits++; }
  }

  if (n_hits == 0) { return false; }

  for (int i = 0; i < e->n_sides(); i++) {
    double temp;
    if (plucker_test(e->side_ptr(i), r, -u, temp)) { n_hits++; }
  }

  // should get at least 2 hits if the
  // location is inside or on a tet boundary
  return n_hits >= 2;
}

void
LibMesh::intersect_track(libMesh::Point start,
                         libMesh::Point dir,
                         double track_len,
                         UnstructuredMeshHits& hits) const
{
  double track_remaining = track_len;

  // attempt to locate a tet for the starting point
  int thread = omp_get_thread_num();
  auto e = (*point_locators_[thread])(start);

  // the point_locator seems off sometimes,
  // ensure the point is actually in the tet
  // and check the neighbors as well
  if (e && !inside_tet(start, dir, e)) {
    bool found = false;
    // try to find a new tet adjacent to this one
    for (int i = 0; i < e->n_neighbors(); i ++) {
      if (e->neighbor_ptr(i) && inside_tet(start, dir, e->neighbor_ptr(i))) {
        e = e->neighbor_ptr(i);
        found = true;
        break;
      }
    }

    if (!found) {
      warning("Starting point is not inside the specified tet.");
    }
  }

  // if the point locator fails, look for mesh
  // entry along the track
  if (!e) {
    auto result = locate_boundary_element(start, start + track_len * dir);
    if (result.second) {
      // if an intersection is found, update the remaining track length
      // and set the element
      e = result.second;
      track_remaining -= result.first;
      // advance position along track
      start += dir * result.first;
    } else {
      // if there was no intersection, we're done
      return;
    }
  }

  std::set<libMesh::dof_id_type> visited;

  bool first = true;

  while (true) {
    // find the positive distance triangle intersection
    double dist = -1.0;
    int side = -1;

    for (int i = 0; i < e->n_sides(); i++) {
      auto tri = e->side_ptr(i);
      if (visited.count(e->key(i))) {
        continue;
      }
      if (tri->type() != libMesh::ElemType::TRI3) { warning("Non-triangle element found"); }
      double temp_dist = -1.0;
      bool hit = plucker_test(e->side_ptr(i), start, dir, temp_dist);
      // if (hit and temp_dist > FP_COINCIDENT) {
      if (hit and temp_dist >= 0 and temp_dist > dist) {
        side = i;
        dist = temp_dist;
        first = false;
      }
    }

    // make sure we found a hit for the tet we're in

    // if we don't find a hit for this tet, we may
    if (side == -1) {

      if (first) {
        warning("Couldn't get hit on first iteration");
        inside_tet(start, dir, e);
        first = false;
      }
      auto orig_e = e;
      start += dir * TINY_BIT; // nudge particle forward

      int thread = omp_get_thread_num();
      e = (*point_locators_[thread])(start);

      if (!e) {
        if (!orig_e->on_boundary()) {
          std::cout << "May have incorrectly truncated a track." << std::endl;
        }
        return;
      }
      continue;
    } else {
      // add hit to output
      hits.push_back(std::pair<double, const libMesh::Elem*>(std::min(track_remaining, dist), e));
      // add this side's centroid to the visited list
      visited.insert(e->key(side));
      // advance position along track
      start += dir * std::min(track_remaining, dist);
      // subtract from remaining track length
      track_remaining -= dist;
    }

    // if we've reached the end of the track, break
    if (track_remaining <= 0.0) { break; }

    // get tet on the other side
    auto next_e = e->neighbor_ptr(side);

    // check to make sure this tet contains our current location
    if (next_e && !next_e->contains_point(start)) {
      warning("Moving into tet that does not contain the current location.");
    }

    // if there is no next element, we may have exited the mesh
    // and will re-enter elsewhere
    if (!next_e) {
      auto result = locate_boundary_element(start,
                                            start + dir * track_remaining);
      if (result.second) {
        // advance position along track
        // to re-entry point
        start += dir * result.first;
        // remove this distance from the remaining track length
        track_remaining -= result.first;
        // update
        e = result.second;
        continue;
      } else {
        // if no next intersection with the mesh, we're done
        return;
      }
    }

    // if the mesh entry point is too far away
    // break
    if (track_remaining <= 0.0) { break; }

    // check that the hit is correct
    if (!elements_share_face(e, next_e, side)) {
      // this should maybe throw an error?
      warning("Incorrect adjacent element found");
      break;
    }

    // update to the element along the track
    e = next_e;
  }
}

bool
LibMesh::elements_share_face(const libMesh::Elem* from,
                             const libMesh::Elem* to,
                             unsigned int side) const
{
  for (auto j : to->side_index_range()) {
    if (from->key(side) == to->key(j)) {
      return true;
    }
  }
  return false;
}

double
LibMesh::first(const libMesh::Node& a,
               const libMesh::Node& b) const
{
  if(a(0) < b(0)) {
    return true;
  } else if(a(0) == b(0)) {
    if(a(1) < b(1)) {
      return true;
    } else if(a(1) == b(1)) {
      if(a(2) < b(2)) {
	return true;
      } else {
        return false;
      }
    } else {
      return false;
    }
  } else {
    return false;
  }
}

void LibMesh::to_hdf5(hid_t group) const
{
  hid_t mesh_group = create_group(group, "mesh " + std::to_string(id_));

  write_dataset(mesh_group, "type", "unstructured");
  write_dataset(mesh_group, "filename", filename_);
  write_dataset(mesh_group, "library", "libmesh");

  // write volume of each tet
  std::vector<double> tet_vols;
  for (int i = 0; i < m_->n_elem(); i++) {
    tet_vols.emplace_back(m_->elem_ref(i).volume());
  }
  write_dataset(mesh_group, "volumes", tet_vols);

  close_group(mesh_group);
}

double
LibMesh::plucker_edge_test(const libMesh::Node& vertexa,
                           const libMesh::Node& vertexb,
                           const libMesh::Point& ray,
                           const libMesh::Point& ray_normal) const
{
  double pip;
  const double near_zero = 10.0 * std::numeric_limits<double>::epsilon();

  if(first(vertexa, vertexb)) {
    const libMesh::Point edge = vertexb - vertexa;
    const libMesh::Point edge_normal = edge.cross(vertexa);
    pip = ray * edge_normal + ray_normal * edge;
  } else {
    const libMesh::Point edge = vertexa-vertexb;
    const libMesh::Point edge_normal = edge.cross(vertexb);
    pip = ray * edge_normal + ray_normal * edge;
    pip = -pip;
  }

  if (near_zero > fabs(pip)) pip = 0.0;

  return pip;
}

/* This test uses the same edge-ray computation for adjacent triangles so that
   rays passing close to edges/nodes are handled consistently.

   Reports intersection type for post processing of special cases. Optionally
   screen by orientation and negative/nonnegative distance limits.

   If screening by orientation, substantial pruning can occur. Indicate
   desired orientation by passing 1 (forward), -1 (reverse), or 0 (no preference).
   Note that triangle orientation is not always the same as surface
   orientation due to non-manifold surfaces.

   N. Platis and T. Theoharis, "Fast Ray-Tetrahedron Intersection using Plücker
   Coordinates", Journal of Graphics Tools, Vol. 8, Part 4, Pages 37-48 (2003). */
bool
LibMesh::plucker_test(std::unique_ptr<const libMesh::Elem> tri,
                      const libMesh::Point& start,
                      const libMesh::Point& dir,
                      double& dist) const
{
  const libMesh::Point raya = dir;
  const libMesh::Point rayb = dir.cross(start);

  // get triangle vertices
  auto node0 = tri->node_ref(0);
  auto node1 = tri->node_ref(1);
  auto node2 = tri->node_ref(2);

  double plucker_coord0 = plucker_edge_test(node0, node1, raya, rayb);
  double plucker_coord1 = plucker_edge_test(node1, node2, raya, rayb);
  if( (0.0<plucker_coord0 && 0.0>plucker_coord1) || (0.0>plucker_coord0 && 0.0<plucker_coord1) ) {
    return false;
  }

  double plucker_coord2 = plucker_edge_test(node2, node0, raya, rayb);
  if( (0.0<plucker_coord1 && 0.0>plucker_coord2) || (0.0>plucker_coord1 && 0.0<plucker_coord2) ||
      (0.0<plucker_coord0 && 0.0>plucker_coord2) || (0.0>plucker_coord0 && 0.0<plucker_coord2) ) {
    return false;
  }

  // check for coplanar case to avoid dividing by zero
  if(0.0==plucker_coord0 && 0.0==plucker_coord1 && 0.0==plucker_coord2) {
    return false;
  }

  // get the distance to intersection
  const double inverse_sum = 1.0/(plucker_coord0+plucker_coord1+plucker_coord2);
  assert(0.0 != inverse_sum);
  const libMesh::Point intersection(plucker_coord0*inverse_sum*node2+
                                    plucker_coord1*inverse_sum*node0+
                                    plucker_coord2*inverse_sum*node1);

  // To minimize numerical error, get index of largest magnitude direction.
  int idx = 0;
  double max_abs_dir = 0;
  for(unsigned int i=0; i<3; ++i) {
    if( fabs(dir(i)) > max_abs_dir ) {
      idx = i;
      max_abs_dir = fabs(dir(i));
    }
  }

  // no negative distances
  double temp_dist = (intersection(idx) - start(idx)) / dir(idx);
  if ( fabs(temp_dist) < TINY_BIT ) { temp_dist = 0.0; }
  if ( temp_dist < 0 ) { return false; }

  dist = temp_dist;

  // dist = (intersection - start).norm();

  return true;
}


#endif // LIBMESH

//==============================================================================
// Non-member functions
//==============================================================================

void read_meshes(pugi::xml_node root)
{
  for (auto node : root.children("mesh")) {
    std::string mesh_type;
    if (check_for_node(node, "type")) {
      mesh_type = get_node_value(node, "type", true, true);
    } else {
      mesh_type = "regular";
    }

    std::string mesh_lib;
    if (check_for_node(node, "library")) {
      mesh_lib = get_node_value(node, "library", true, true);
    } else {
      mesh_lib = "moab";
    }

    // Read mesh and add to vector
    if (mesh_type == "regular") {
      model::meshes.push_back(std::make_unique<RegularMesh>(node));
    } else if (mesh_type == "rectilinear") {
      model::meshes.push_back(std::make_unique<RectilinearMesh>(node));
#ifdef DAGMC
    }
    else if (mesh_type == "unstructured" && mesh_lib == "moab") {
      model::meshes.push_back(std::make_unique<UnstructuredMesh>(node));
#endif
#ifdef LIBMESH
    }
    else if (mesh_type == "unstructured" && mesh_lib == "libmesh") {
      std::cout << "Making libmesh mesh" << std::endl;
      model::meshes.push_back(std::make_unique<LibMesh>(node));
#endif
    } else {
      fatal_error("Invalid mesh type: " + mesh_type);
    }

    // Map ID to position in vector
    model::mesh_map[model::meshes.back()->id_] = model::meshes.size() - 1;
  }
}

void meshes_to_hdf5(hid_t group)
{
  // Write number of meshes
  hid_t meshes_group = create_group(group, "meshes");
  int32_t n_meshes = model::meshes.size();
  write_attribute(meshes_group, "n_meshes", n_meshes);

  if (n_meshes > 0) {
    // Write IDs of meshes
    std::vector<int> ids;
    for (const auto& m : model::meshes) {
      m->to_hdf5(meshes_group);
      ids.push_back(m->id_);
    }
    write_attribute(meshes_group, "ids", ids);
  }

  close_group(meshes_group);
}

void free_memory_mesh()
{
  model::meshes.clear();
  model::mesh_map.clear();
}

extern "C" int n_meshes() { return model::meshes.size(); }

} // namespace openmc
