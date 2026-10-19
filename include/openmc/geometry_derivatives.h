#ifndef OPENMC_GEOMETRY_DERIVATIVES_H
#define OPENMC_GEOMETRY_DERIVATIVES_H

#include <cstddef>
#include <cstdint>
#include <unordered_map>
#include <unordered_set>

#include "openmc/memory.h"
#include "openmc/particle.h"
#include "openmc/tensor.h"
#include "openmc/vector.h"
#include "openmc/xml_interface.h"
namespace openmc {

class GeometryDerivative;

namespace model {
extern std::unordered_map<int32_t, int32_t> geometry_derivatives_map;
extern vector<unique_ptr<GeometryDerivative>> geometry_derivatives;
extern std::unordered_set<int32_t> geometry_derivative_tallies;
} // namespace model

void read_geometry_derivatives(pugi::xml_node node);
void prepare_geometry_derivatives();
void tally_geometry_derivatives(
  Particle& p, int32_t tally_id, int64_t score_bin, double score);
void accumulate_geometry_derivatives();
void report_geometry_derivatives();
void reset_geometry_derivatives();
void free_memory_geometry_derivatives();

class GeometryDerivative {

public:
  // Constructors
  GeometryDerivative(pugi::xml_node node);

  // Methods
  void init_results();
  void reset();

  // Accessors
  int32_t id() const { return id_; }
  int32_t tally_id() const { return tally_id_; }
  const auto& surface_ids() const { return surface_ids_; }
  const auto& surface_parameters() const { return surface_parameters_; }

  auto& tally_derivatives() { return tally_derivatives_; }
  const auto& tally_derivatives() const { return tally_derivatives_; }
  const auto& surface_parameter_map() const { return surface_parameter_map_; }
  void accumulate();

private:
  int32_t id_;
  int32_t tally_id_;
  vector<int32_t> surface_ids_;
  vector<std::string> surface_parameters_;


  //! Type alias and custom hash function for surface-parameter pair
  using SurfaceParameterPair = std::pair<int, int>;
  struct SurfaceParameterPairHash {
    std::size_t operator()(const SurfaceParameterPair& pair) const {
      return std::hash<int>()(pair.first) ^ std::hash<int>()(pair.second);
    }
  };

  //! Map surface IDs and parameter indices to flattened result indices.
  //! The first key is the surface ID and the second key is the parameter index.
  std::unordered_map<SurfaceParameterPair, int32_t, SurfaceParameterPairHash> surface_parameter_map_;

  //! Results of the tally derivative with respsect to various geometric
  //! parameters. The first dimension of the array is for geometric parameter
  //! indices. The second dimension is the flattened tally bin index. The
  //! third dimension is the tally for the current batch of the tally derivative
  //! with respect to the corresponding tally bin and geometric parameter for
  //! the current batch. The fourth dimension is the final result with
  //! accumulation after each batch.
  tensor::Tensor<double> tally_derivatives_;
};

} // namespace openmc

extern "C" {

int openmc_get_geometry_derivative_index(int32_t id, int32_t* index);
int openmc_geometry_derivative_get_id(int32_t index, int32_t* id);
int openmc_geometry_derivative_get_tally_id(int32_t index, int32_t* id);
int openmc_geometry_derivative_get_surface_ids(
  int32_t index, const int32_t** ids, size_t* n);
int openmc_geometry_derivative_get_surface_parameters(
  int32_t index, int32_t surface_id, const char** parameters);
int openmc_geometry_derivative_results(
  int32_t index, double** results, size_t* shape);
int openmc_geometry_derivative_reset(int32_t index);
size_t openmc_geometry_derivatives_size();
}

#endif // OPENMC_GEOMETRY_DERIVATIVES_H
