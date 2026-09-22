#include "openmc/geometry_derivatives.h"

#include "openmc/cell.h"
#include "openmc/error.h"
#include "openmc/container_util.h"

#include "openmc/tallies/tally.h"

namespace openmc {

namespace model {
  std::unordered_map<int32_t, int32_t> geometry_derivatives_map;
  vector<unique_ptr<GeometryDerivative>> geometry_derivatives;
  std::unordered_set<int32_t> derivative_surface_indices;
}

GeometryDerivative::GeometryDerivative(pugi::xml_node node)
{
  if (!check_for_node(node, "id")) {
    fatal_error("Must specify id of geometry derivative in geometry_derivatives XML file.");
  }
  id_ = std::stoi(get_node_value(node, "id"));

  if (!check_for_node(node, "tally")) {
    fatal_error("Must specify tally_id of geometry derivative in geometry_derivatives XML file.");
  }
  tally_id_ = std::stoi(get_node_value(node, "tally"));

  if (!check_for_node(node, "cell")) {
    fatal_error("Must specify cell of geometry derivative in geometry_derivatives XML file.");
  }
  cell_id_ = std::stoi(get_node_value(node, "cell"));
}

void GeometryDerivative::init_results()
{
  surface_indices_.clear();
  geom_parameters_ = {};
  tally_derivatives_ = {};

  // determine number of differentiable parameters on the cell
  const auto& cell = model::cells[model::cell_map[cell_id_]];
  const auto& surfaces = cell->surfaces();

  int n_surface_params = 0;
  for (const auto& surf_token : surfaces) {
    int32_t surface_index = std::abs(surf_token) - 1;
    model::derivative_surface_indices.insert(surface_index);
    const auto& surf = model::surfaces[surface_index];
    surface_indices_[surf->id_] = n_surface_params;
    n_surface_params += surf->n_diff_params();
  }

  const auto& tally = model::tallies[model::tally_map[tally_id_]];
  int n_scores = tally->n_scores();

  geom_parameters_ = tensor::Tensor<double>(
    {static_cast<size_t>(n_surface_params), 4});

  tally_derivatives_ = tensor::Tensor<double>(
    {static_cast<size_t>(n_surface_params), static_cast<size_t>(n_scores)});
}

void read_geometry_derivatives(pugi::xml_node node)
{
  if (!check_for_node(node, "geometry_derivative")) {
    return;
  }

  for (pugi::xml_node geom_deriv_node : node.children("geometry_derivative")) {
    model::geometry_derivatives.push_back(std::make_unique<GeometryDerivative>(geom_deriv_node));
    model::geometry_derivatives_map[model::geometry_derivatives.back()->id()] = model::geometry_derivatives.size() - 1;
  }
}

void prepare_geometry_derivatives()
{
  model::derivative_surface_indices.clear();
  for (const auto& geom_deriv : model::geometry_derivatives) {
    geom_deriv->init_results();
  }
}

void update_surface_derivative(Particle& p)
{
  // compute the derivatives for the surface parameters
  // of the surface being crossed
  const auto& surface = model::surfaces[p.boundary().surface_index()];
  std::vector<double> surface_derivatives = surface->derivatives(p.r(), p.u());
  int32_t surface_id = surface->id_;
  // find the geometry derivative for the cell being crossed
  for (const auto& geom_deriv : model::geometry_derivatives) {
    // find the index of the surface in the geometry derivative
    if (geom_deriv->surface_indices().count(surface_id) > 0) {
      int32_t surface_param_start = geom_deriv->surface_indices()[surface_id];
      for (int i = 0; i < surface_derivatives.size(); ++i) {
        // update jacobian and jacobian derivative for the surface parameter
#pragma omp atomic
        geom_deriv->geom_parameters()(surface_param_start + i, 0) *= p.boundary().distance();
#pragma omp atomic
        geom_deriv->geom_parameters()(surface_param_start + i, 1) += surface_derivatives[i]  / p.boundary().distance();
        // update attenuation derivative factor for the surface parameter
#pragma omp atomic
        geom_deriv->geom_parameters()(surface_param_start + i, 2) *= std::exp(-p.macro_xs().total * p.collision_distance());
        double new_att_deriv = -p.macro_xs().total * p.collision_distance() / p.boundary().distance() * surface_derivatives[i];
#pragma omp atomic
        geom_deriv->geom_parameters()(surface_param_start + i, 3) += new_att_deriv;
      }
    }
  }
}

} // namespace openmc