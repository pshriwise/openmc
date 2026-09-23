#include "openmc/geometry_derivatives.h"

#include "openmc/capi.h"
#include "openmc/cell.h"
#include "openmc/container_util.h"
#include "openmc/error.h"
#include "openmc/settings.h"
#include "openmc/simulation.h"
#include "openmc/source.h"

#include "openmc/tallies/tally.h"

namespace openmc {

namespace model {
std::unordered_map<int32_t, int32_t> geometry_derivatives_map;
vector<unique_ptr<GeometryDerivative>> geometry_derivatives;
std::unordered_set<int32_t> derivative_surface_indices;
std::unordered_set<int32_t> geometry_derivative_tallies;
} // namespace model

GeometryDerivative::GeometryDerivative(pugi::xml_node node)
{
  if (!check_for_node(node, "id")) {
    fatal_error("Must specify id of geometry derivative in "
                "geometry_derivatives XML file.");
  }
  id_ = std::stoi(get_node_value(node, "id"));

  if (!check_for_node(node, "tally")) {
    fatal_error("Must specify tally_id of geometry derivative in "
                "geometry_derivatives XML file.");
  }
  tally_id_ = std::stoi(get_node_value(node, "tally"));

  if (!check_for_node(node, "cell")) {
    fatal_error("Must specify cell of geometry derivative in "
                "geometry_derivatives XML file.");
  }
  cell_id_ = std::stoi(get_node_value(node, "cell"));
}

void GeometryDerivative::init_results()
{
  surface_indices_.clear();
  surface_ids_.clear();
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
    surface_ids_.push_back(surf->id_);
    n_surface_params += surf->n_diff_params();
  }

  const auto& tally = model::tallies[model::tally_map[tally_id_]];

  // ensure that the tally is using a collision estimator (for now)
  if (tally->estimator_ != TallyEstimator::COLLISION) {
    fatal_error(fmt::format(
      "Tally {} is not using a collision estimator. Geometry derivatives "
      "are only supported for collision estimators.",
      tally_id_));
  }

  geom_parameters_ =
    tensor::Tensor<double>({static_cast<size_t>(n_surface_params), 4});

  int n_tally_bins = tally->results().shape(0) * tally->results().shape(1);
  tally_derivatives_ =
    tensor::Tensor<double>({static_cast<size_t>(n_surface_params),
      static_cast<size_t>(n_tally_bins), 2});
}

void GeometryDerivative::reset()
{
  if (geom_parameters_.size() != 0) {
    geom_parameters_.fill(0.0);
  }
  if (tally_derivatives_.size() != 0) {
    tally_derivatives_.fill(0.0);
  }
}

void GeometryDerivative::accumulate()
{
  double total_source = settings::run_mode == RunMode::FIXED_SOURCE
                          ? model::external_sources_probability.integral()
                          : 1.0;
  double contributing_particles = settings::reduce_tallies
                                    ? settings::n_particles
                                    : simulation::work_per_rank;
  double norm =
    total_source / (contributing_particles * settings::gen_per_batch);
  if (settings::solver_type == SolverType::RANDOM_RAY) {
    norm = 1.0;
  }

  for (int i = 0; i < geom_parameters_.shape(0); ++i) {
    for (int j = 0; j < tally_derivatives_.shape(1); ++j) {
      tally_derivatives_(i, j, 1) += tally_derivatives_(i, j, 0) * norm;
      tally_derivatives_(i, j, 0) = 0.0;
    }
  }
}

void read_geometry_derivatives(pugi::xml_node node)
{
  if (!check_for_node(node, "geometry_derivative")) {
    return;
  }

  for (pugi::xml_node geom_deriv_node : node.children("geometry_derivative")) {
    model::geometry_derivatives.push_back(
      std::make_unique<GeometryDerivative>(geom_deriv_node));
    model::geometry_derivatives_map[model::geometry_derivatives.back()->id()] =
      model::geometry_derivatives.size() - 1;
  }
}

void prepare_geometry_derivatives()
{
  model::derivative_surface_indices.clear();
  model::geometry_derivative_tallies.clear();
  for (const auto& geom_deriv : model::geometry_derivatives) {
    geom_deriv->init_results();
    model::geometry_derivative_tallies.insert(geom_deriv->tally_id());
  }
}

void update_surface_derivative(Particle& p)
{
  // compute the derivatives for the surface parameters
  // of the surface being crossed
  const auto& surface = model::surfaces[p.boundary().surface_index()];
  std::vector<double> surface_derivatives = surface->derivatives(p.r(), p.u());
  int32_t surface_id = surface->id_;
  double adv_distance =
    std::min(p.boundary().distance(), p.collision_distance());
  // find the geometry derivative for the cell being crossed
  for (const auto& geom_deriv : model::geometry_derivatives) {
    // find the index of the surface in the geometry derivative
    if (geom_deriv->surface_indices().count(surface_id) > 0) {
      int32_t surface_param_start = geom_deriv->surface_indices()[surface_id];
      for (int i = 0; i < surface_derivatives.size(); ++i) {
        // update attenuation derivative factor for the surface parameter
#pragma omp atomic
        geom_deriv->geom_parameters()(surface_param_start + i, 2) *=
          std::exp(-p.macro_xs().total * adv_distance);
        double new_att_deriv = -p.macro_xs().total * adv_distance /
                               p.boundary().distance() * surface_derivatives[i];
#pragma omp atomic
        geom_deriv->geom_parameters()(surface_param_start + i, 3) +=
          new_att_deriv;
        if (p.boundary().distance() <= p.collision_distance()) {
          continue;
        }
        // update jacobian and jacobian derivative for the surface parameter
#pragma omp atomic
        geom_deriv->geom_parameters()(surface_param_start + i, 0) *=
          p.boundary().distance();
#pragma omp atomic
        geom_deriv->geom_parameters()(surface_param_start + i, 1) +=
          surface_derivatives[i] / p.boundary().distance();
      }
    }
  }
}

void tally_geometry_derivatives(
  int32_t tally_id, int64_t score_bin, double score)
{
  for (const auto& geom_deriv : model::geometry_derivatives) {
    if (geom_deriv->tally_id() != tally_id) {
      continue;
    }

    for (int i = 0; i < geom_deriv->geom_parameters().shape(0); ++i) {
      // update the tally derivative with respect to the surface parameter
      double dJ = geom_deriv->geom_parameters()(i, 1);
      double datt = geom_deriv->geom_parameters()(i, 3);
      double deriv = (dJ + datt) * score;
#pragma omp atomic
      geom_deriv->tally_derivatives()(i, score_bin, 0) += deriv;
    }
  }
}

void accumulate_geometry_derivatives()
{
  for (const auto& geom_deriv : model::geometry_derivatives) {
    geom_deriv->accumulate();
  }
}

void report_geometry_derivatives()
{
  return;
  for (const auto& geom_deriv : model::geometry_derivatives) {
    std::cout << "Geometry Derivative ID: " << geom_deriv->id() << std::endl;
    std::cout << "Tally ID: " << geom_deriv->tally_id() << std::endl;
    std::cout << "Cell ID: " << geom_deriv->cell_id() << std::endl;
    std::cout << "Surface Indices: ";
    for (const auto& [surface_id, surface_param_start] :
      geom_deriv->surface_indices()) {
      std::cout << surface_id << " ";
    }
    std::cout << std::endl;
    // report geometry parameters and tally derivatives in loops
    for (int i = 0; i < geom_deriv->geom_parameters().shape(0); ++i) {
      std::cout << "Surface Parameter Index: " << i << std::endl;
      std::cout << "Jacobian: " << geom_deriv->geom_parameters()(i, 0)
                << std::endl;
      std::cout << "Jacobian Derivative: "
                << geom_deriv->geom_parameters()(i, 1) << std::endl;
      std::cout << "Attenuation Factor: " << geom_deriv->geom_parameters()(i, 2)
                << std::endl;
      std::cout << "Attenuation Factor Derivative: "
                << geom_deriv->geom_parameters()(i, 3) << std::endl;
    }
    for (int i = 0; i < geom_deriv->tally_derivatives().shape(0); ++i) {
      std::cout << "Tally Derivative for Surface Parameter Index: " << i
                << std::endl;
      for (int j = 0; j < geom_deriv->tally_derivatives().shape(1); ++j) {
        std::cout << "Score Bin " << j << ": "
                  << geom_deriv->tally_derivatives()(i, j, 1) << std::endl;
      }
    }
    std::cout << "----------------------------------------" << std::endl;
  }
}

void reset_geometry_derivatives()
{
  for (const auto& geom_deriv : model::geometry_derivatives) {
    geom_deriv->reset();
  }
}

void free_memory_geometry_derivatives()
{
  model::geometry_derivatives_map.clear();
  model::geometry_derivatives.clear();
  model::derivative_surface_indices.clear();
  model::geometry_derivative_tallies.clear();
}

namespace {

GeometryDerivative* get_geometry_derivative(int32_t index)
{
  if (index < 0 || index >= model::geometry_derivatives.size()) {
    set_errmsg("Index in geometry derivatives array is out of bounds.");
    return nullptr;
  }
  return model::geometry_derivatives[index].get();
}

} // namespace

} // namespace openmc

using namespace openmc;

extern "C" int openmc_get_geometry_derivative_index(int32_t id, int32_t* index)
{
  auto it = model::geometry_derivatives_map.find(id);
  if (it == model::geometry_derivatives_map.end()) {
    set_errmsg(fmt::format("No geometry derivative exists with ID={}.", id));
    return OPENMC_E_INVALID_ID;
  }

  *index = it->second;
  return 0;
}

extern "C" int openmc_geometry_derivative_get_id(int32_t index, int32_t* id)
{
  auto* deriv = get_geometry_derivative(index);
  if (!deriv)
    return OPENMC_E_OUT_OF_BOUNDS;

  *id = deriv->id();
  return 0;
}

extern "C" int openmc_geometry_derivative_get_tally_id(
  int32_t index, int32_t* id)
{
  auto* deriv = get_geometry_derivative(index);
  if (!deriv)
    return OPENMC_E_OUT_OF_BOUNDS;

  *id = deriv->tally_id();
  return 0;
}

extern "C" int openmc_geometry_derivative_get_cell_id(
  int32_t index, int32_t* id)
{
  auto* deriv = get_geometry_derivative(index);
  if (!deriv)
    return OPENMC_E_OUT_OF_BOUNDS;

  *id = deriv->cell_id();
  return 0;
}

extern "C" int openmc_geometry_derivative_get_surface_ids(
  int32_t index, const int32_t** surface_ids, size_t* n)
{
  auto* deriv = get_geometry_derivative(index);
  if (!deriv)
    return OPENMC_E_OUT_OF_BOUNDS;

  *surface_ids = deriv->surface_ids().data();
  *n = deriv->surface_ids().size();
  return 0;
}

extern "C" int openmc_geometry_derivative_parameters(
  int32_t index, double** parameters, size_t* shape)
{
  auto* deriv = get_geometry_derivative(index);
  if (!deriv)
    return OPENMC_E_OUT_OF_BOUNDS;

  auto& geom_parameters = deriv->geom_parameters();
  if (geom_parameters.size() == 0) {
    set_errmsg("Geometry derivative parameters have not been allocated yet.");
    return OPENMC_E_ALLOCATE;
  }

  *parameters = geom_parameters.data();
  auto s = geom_parameters.shape();
  shape[0] = s[0];
  shape[1] = s[1];
  return 0;
}

extern "C" int openmc_geometry_derivative_results(
  int32_t index, double** results, size_t* shape)
{
  auto* deriv = get_geometry_derivative(index);
  if (!deriv)
    return OPENMC_E_OUT_OF_BOUNDS;

  auto& tally_derivatives = deriv->tally_derivatives();
  if (tally_derivatives.size() == 0) {
    set_errmsg("Geometry derivative results have not been allocated yet.");
    return OPENMC_E_ALLOCATE;
  }

  *results = tally_derivatives.data();
  auto s = tally_derivatives.shape();
  shape[0] = s[0];
  shape[1] = s[1];
  shape[2] = s[2];
  return 0;
}

extern "C" int openmc_geometry_derivative_reset(int32_t index)
{
  auto* deriv = get_geometry_derivative(index);
  if (!deriv)
    return OPENMC_E_OUT_OF_BOUNDS;

  deriv->reset();
  return 0;
}

extern "C" size_t openmc_geometry_derivatives_size()
{
  return model::geometry_derivatives.size();
}
