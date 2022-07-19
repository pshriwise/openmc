#include "openmc/distribution_spatial.h"

#include "openmc/error.h"
#include "openmc/random_lcg.h"
#include "openmc/xml_interface.h"
#include "openmc/mesh.h"
#include "openmc/search.h"

#include <iostream>
#include <fstream>

namespace openmc {

//==============================================================================
// CartesianIndependent implementation
//==============================================================================

CartesianIndependent::CartesianIndependent(pugi::xml_node node)
{
  // Read distribution for x coordinate
  if (check_for_node(node, "x")) {
    pugi::xml_node node_dist = node.child("x");
    x_ = distribution_from_xml(node_dist);
  } else {
    // If no distribution was specified, default to a single point at x=0
    double x[] {0.0};
    double p[] {1.0};
    x_ = UPtrDist {new Discrete {x, p, 1}};
  }

  // Read distribution for y coordinate
  if (check_for_node(node, "y")) {
    pugi::xml_node node_dist = node.child("y");
    y_ = distribution_from_xml(node_dist);
  } else {
    // If no distribution was specified, default to a single point at y=0
    double x[] {0.0};
    double p[] {1.0};
    y_ = UPtrDist {new Discrete {x, p, 1}};
  }

  // Read distribution for z coordinate
  if (check_for_node(node, "z")) {
    pugi::xml_node node_dist = node.child("z");
    z_ = distribution_from_xml(node_dist);
  } else {
    // If no distribution was specified, default to a single point at z=0
    double x[] {0.0};
    double p[] {1.0};
    z_ = UPtrDist {new Discrete {x, p, 1}};
  }
}

Position CartesianIndependent::sample(uint64_t* seed) const
{
  return {x_->sample(seed), y_->sample(seed), z_->sample(seed)};
}

//==============================================================================
// CylindricalIndependent implementation
//==============================================================================

CylindricalIndependent::CylindricalIndependent(pugi::xml_node node)
{
  // Read distribution for r-coordinate
  if (check_for_node(node, "r")) {
    pugi::xml_node node_dist = node.child("r");
    r_ = distribution_from_xml(node_dist);
  } else {
    // If no distribution was specified, default to a single point at r=0
    double x[] {0.0};
    double p[] {1.0};
    r_ = make_unique<Discrete>(x, p, 1);
  }

  // Read distribution for phi-coordinate
  if (check_for_node(node, "phi")) {
    pugi::xml_node node_dist = node.child("phi");
    phi_ = distribution_from_xml(node_dist);
  } else {
    // If no distribution was specified, default to a single point at phi=0
    double x[] {0.0};
    double p[] {1.0};
    phi_ = make_unique<Discrete>(x, p, 1);
  }

  // Read distribution for z-coordinate
  if (check_for_node(node, "z")) {
    pugi::xml_node node_dist = node.child("z");
    z_ = distribution_from_xml(node_dist);
  } else {
    // If no distribution was specified, default to a single point at z=0
    double x[] {0.0};
    double p[] {1.0};
    z_ = make_unique<Discrete>(x, p, 1);
  }

  // Read cylinder center coordinates
  if (check_for_node(node, "origin")) {
    auto origin = get_node_array<double>(node, "origin");
    if (origin.size() == 3) {
      origin_ = origin;
    } else {
      fatal_error(
        "Origin for cylindrical source distribution must be length 3");
    }
  } else {
    // If no coordinates were specified, default to (0, 0, 0)
    origin_ = {0.0, 0.0, 0.0};
  }
}

Position CylindricalIndependent::sample(uint64_t* seed) const
{
  double r = r_->sample(seed);
  double phi = phi_->sample(seed);
  double x = r * cos(phi) + origin_.x;
  double y = r * sin(phi) + origin_.y;
  double z = z_->sample(seed) + origin_.z;
  return {x, y, z};
}

//==============================================================================
// SphericalIndependent implementation
//==============================================================================

SphericalIndependent::SphericalIndependent(pugi::xml_node node)
{
  // Read distribution for r-coordinate
  if (check_for_node(node, "r")) {
    pugi::xml_node node_dist = node.child("r");
    r_ = distribution_from_xml(node_dist);
  } else {
    // If no distribution was specified, default to a single point at r=0
    double x[] {0.0};
    double p[] {1.0};
    r_ = make_unique<Discrete>(x, p, 1);
  }

  // Read distribution for cos_theta-coordinate
  if (check_for_node(node, "cos_theta")) {
    pugi::xml_node node_dist = node.child("cos_theta");
    cos_theta_ = distribution_from_xml(node_dist);
  } else {
    // If no distribution was specified, default to a single point at cos_theta=0
    double x[] {0.0};
    double p[] {1.0};
    cos_theta_ = make_unique<Discrete>(x, p, 1);
  }

  // Read distribution for phi-coordinate
  if (check_for_node(node, "phi")) {
    pugi::xml_node node_dist = node.child("phi");
    phi_ = distribution_from_xml(node_dist);
  } else {
    // If no distribution was specified, default to a single point at phi=0
    double x[] {0.0};
    double p[] {1.0};
    phi_ = make_unique<Discrete>(x, p, 1);
  }

  // Read sphere center coordinates
  if (check_for_node(node, "origin")) {
    auto origin = get_node_array<double>(node, "origin");
    if (origin.size() == 3) {
      origin_ = origin;
    } else {
      fatal_error("Origin for spherical source distribution must be length 3");
    }
  } else {
    // If no coordinates were specified, default to (0, 0, 0)
    origin_ = {0.0, 0.0, 0.0};
  }
}

Position SphericalIndependent::sample(uint64_t* seed) const
{
  double r = r_->sample(seed);
  double cos_theta = cos_theta_->sample(seed);
  double phi = phi_->sample(seed);
  // sin(theta) by sin**2 + cos**2 = 1
  double x = r * std::sqrt(1 - cos_theta * cos_theta) * cos(phi) + origin_.x;
  double y = r * std::sqrt(1 - cos_theta * cos_theta) * sin(phi) + origin_.y;
  double z = r * cos_theta + origin_.z;
  return {x, y, z};
}

//==============================================================================
// MeshIndependent implementation
//==============================================================================

// Loosely adapted Patrick's code shown here
MeshIndependent::MeshIndependent(pugi::xml_node node)
{
  // No in-tet distributions implemented, could include distributions for the barycentric coords
  // Read in unstructured mesh from mesh_id value
  int32_t mesh_id = std::stoi(get_node_value(node, "mesh_id"));
  mesh_map_idx_ = model::mesh_map.at(mesh_id);
  const auto& mesh_ptr = model::meshes[mesh_map_idx_];
  
  // Check whether mesh pointer points to an unstructured mesh
  umesh_ptr_ = dynamic_cast<UnstructuredMesh*>(mesh_ptr.get());
  if (!umesh_ptr_) {fatal_error("Mesh passed to spatial distribution is not an unstructured mesh object"); }

  // Create CDF based on weighting scheme
  tot_bins_ = umesh_ptr_->n_bins();
  // float weights[tot_bins_] = {};
  std::vector<double> weights = {};
  weights.resize(tot_bins_);
  double temp_total_weight = 0.0;
  mesh_CDF_.resize(tot_bins_+1);
  mesh_CDF_[0] = {0.0};
  total_weight_ = 0.0;
  mesh_weights_.resize(get_node_array<double>(node, "weights_from_file").size());

  // Create cdfs for sampling for an element over a mesh
  // Volume scheme is weighted by the volume of each tet
  // File scheme is weighted by an array given in the xml file

  if (check_for_node(node, "elem_weight_scheme")) {
    sample_scheme_ = get_node_value(node, "elem_weight_scheme");
  } else {
      fatal_error("No weighting scheme was specified for element sampling over the mesh");
  }

  if (sample_scheme_ == "volume"){
    for (int i = 0; i<tot_bins_; i++){
      weights[i] = umesh_ptr_->volume(i);
    }
  } else if (sample_scheme_ == "file"){

    if (check_for_node(node, "weights_from_file")) {
      mesh_weights_ = get_node_array<double>(node, "weights_from_file");
      if (mesh_weights_.size() != tot_bins_){
        write_message("The size of the weights array from the xml file does not equal the number of elements in the mesh.");
        // fatal_error("The size of the weights array from the xml file does not equal the number of elements in the mesh.");
      }
      for (int i = 0; i < tot_bins_; i++){
        weights[i] = mesh_weights_[i];
      } 
    } else {
        fatal_error("No weights given in xml file");
    }
  } else{
    fatal_error("Type of element weighting scheme provided is not in the supported types (volume or file)");
  }
  for (int i = 0; i<tot_bins_; i++){
    temp_total_weight = temp_total_weight + weights[i];
  }
  total_weight_ = temp_total_weight;
  for (int i = 0; i<tot_bins_; i++){
    mesh_CDF_[i+1] = mesh_CDF_[i] + weights[i]/total_weight_;
  }
}

Position MeshIndependent::sample(uint64_t* seed) const
{ 
  // Create random variable for sampling element from mesh
  float eta = prn(seed);
  // Sample over the CDF defined in initialization above
  int32_t tet_bin = lower_bound_index(mesh_CDF_.begin(), mesh_CDF_.end(), eta);
  return umesh_ptr_->sample(seed, tet_bin);
}


//==============================================================================
// SpatialBox implementation
//==============================================================================

SpatialBox::SpatialBox(pugi::xml_node node, bool fission)
  : only_fissionable_ {fission}
{
  // Read lower-right/upper-left coordinates
  auto params = get_node_array<double>(node, "parameters");
  if (params.size() != 6)
    openmc::fatal_error("Box/fission spatial source must have six "
                        "parameters specified.");

  lower_left_ = Position {params[0], params[1], params[2]};
  upper_right_ = Position {params[3], params[4], params[5]};
}

Position SpatialBox::sample(uint64_t* seed) const
{
  Position xi {prn(seed), prn(seed), prn(seed)};
  return lower_left_ + xi * (upper_right_ - lower_left_);
}

//==============================================================================
// SpatialPoint implementation
//==============================================================================

SpatialPoint::SpatialPoint(pugi::xml_node node)
{
  // Read location of point source
  auto params = get_node_array<double>(node, "parameters");
  if (params.size() != 3)
    openmc::fatal_error("Point spatial source must have three "
                        "parameters specified.");

  // Set position
  r_ = Position {params.data()};
}

Position SpatialPoint::sample(uint64_t* seed) const
{
  return r_;
}

} // namespace openmc
