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
  extern std::unordered_set<int32_t> derivative_surface_indices;
}

void read_geometry_derivatives(pugi::xml_node node);
void prepare_geometry_derivatives();
void update_surface_derivative(Particle& p);

class GeometryDerivative {

  public:
  // Constructors
  GeometryDerivative(pugi::xml_node node);

  // Methods
  void init_results();

  // Accessors
  int32_t id() const { return id_; }
  int32_t tally_id() const { return tally_id_; }
  int32_t cell_id() const { return cell_id_; }

  const auto& surface_indices() const { return surface_indices_; }
  auto& surface_indices() { return surface_indices_; }
  const auto& goem_parameters() const { return geom_parameters_; }
  auto& geom_parameters() { return geom_parameters_; }
  auto& tally_derivatives() { return tally_derivatives_; }
  const auto& tally_derivatives() const { return tally_derivatives_; }


  private:
  int32_t id_;
  int32_t tally_id_;
  int32_t cell_id_;

  //! Results of the geometry derivative tally -- the first dimesion of the
  //! array is for the surface index. The second dimension is size 5 holding the
  //! following values for each geometric paramater/tally score combination:
  //! 0: jacobian
  //! 1: jacobian derivative
  //! 2: attenuation factor
  //! 3: attenuation factor derivative
  tensor::Tensor<double> geom_parameters_;

  //! Results of the tally derivative with respsect to various geometric
  //! parameters. The first dimension of the array is for geometric parameter
  //! indices. The second dimension is the flattened tally bin index. The
  //! third dimension is the result of the tally derivative with respect to
  //! the corresponding tally bin and geometric parameter.
  tensor::Tensor<double> tally_derivatives_;

  //! Store start in surface parameter index for the surface of each cell
  std::unordered_map<int32_t, int32_t> surface_indices_;
};

} // namespace openmc