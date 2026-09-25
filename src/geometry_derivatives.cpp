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

  int n_tally_bins = tally->results().shape(0) * tally->results().shape(1);
  tally_derivatives_ =
    tensor::Tensor<double>({static_cast<size_t>(n_surface_params),
      static_cast<size_t>(n_tally_bins), 6});
}

void GeometryDerivative::reset()
{
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

  for (int i = 0; i < tally_derivatives_.shape(0); ++i) {
    for (int j = 0; j < tally_derivatives_.shape(1); ++j) {
      tally_derivatives_(i, j, 1) += tally_derivatives_(i, j, 0) * norm;
      tally_derivatives_(i, j, 0) = 0.0;
      tally_derivatives_(i, j, 3) += tally_derivatives_(i, j, 2) * norm;
      tally_derivatives_(i, j, 2) = 0.0;
      tally_derivatives_(i, j, 5) += tally_derivatives_(i, j, 4) * norm;
      tally_derivatives_(i, j, 4) = 0.0;
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

void tally_geometry_derivative(int32_t geometry_deriv_idx, int64_t score_bin, int32_t param_idx, double score, double dJ, double datt) {
  const auto& geom_deriv = model::geometry_derivatives[geometry_deriv_idx];
  // update the tally derivative with respect to the surface parameter
  double deriv = (dJ + datt) * score;
#pragma omp atomic
    geom_deriv->tally_derivatives()(param_idx, score_bin, 0) += deriv;
#pragma omp atomic
    geom_deriv->tally_derivatives()(param_idx, score_bin, 2) += score * datt;
#pragma omp atomic
    geom_deriv->tally_derivatives()(param_idx, score_bin, 4) += score * dJ;
}

void tally_geometry_derivatives(Particle& p, int32_t tally_id, int64_t score_bin, double score)
{
  for (const auto& geom_deriv_state : p.geometry_deriv_state()) {
    const auto& geom_deriv = model::geometry_derivatives[geom_deriv_state.geometry_derivative_idx];
    if (geom_deriv->tally_id() != tally_id) {
      continue;
    }
    tally_geometry_derivative(geom_deriv_state.geometry_derivative_idx, score_bin, geom_deriv_state.parameter_idx, score, geom_deriv_state.dj, geom_deriv_state.df);
  }
}

void accumulate_geometry_derivatives()
{
  for (const auto& geom_deriv : model::geometry_derivatives) {
    geom_deriv->accumulate();
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
