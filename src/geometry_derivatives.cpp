#include "openmc/geometry_derivatives.h"

#include "openmc/capi.h"
#include "openmc/cell.h"
#include "openmc/container_util.h"
#include "openmc/error.h"
#include "openmc/settings.h"
#include "openmc/simulation.h"
#include "openmc/source.h"
#include "openmc/string_utils.h"

#include "openmc/tallies/tally.h"

namespace openmc {

namespace model {
std::unordered_map<int32_t, int32_t> geometry_derivatives_map;
vector<unique_ptr<GeometryDerivative>> geometry_derivatives;
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

  for (auto surface_node : node.children("surface")) {
    if (!check_for_node(surface_node, "id") ||
        !check_for_node(surface_node, "parameters")) {
      fatal_error("Must specify id and parameters for each geometry derivative "
                  "surface.");
    }
    int32_t surface_id = std::stoi(get_node_value(surface_node, "id"));
    if (contains(surface_ids_, surface_id)) {
      fatal_error(fmt::format(
        "Duplicate surface {} in geometry derivative {}.", surface_id, id_));
    }
    surface_ids_.push_back(surface_id);
    surface_parameters_.push_back(get_node_value(surface_node, "parameters"));
  }
}

void GeometryDerivative::init_results()
{
  surface_parameter_map_.clear();
  tally_derivatives_ = {};

  int n_params = 0;
  for (int i = 0; i < surface_ids_.size(); ++i) {
    int32_t surface_id = surface_ids_[i];
    auto surface_it = model::surface_map.find(surface_id);
    if (surface_it == model::surface_map.end()) {
      fatal_error(
        fmt::format("Could not find surface {} for geometry derivative {}.",
          surface_id, id_));
    }
    const auto& surface = model::surfaces[surface_it->second];
    for (const auto& param : split(surface_parameters_[i])) {
      int param_idx = surface->diff_param_index(param);
      if (param_idx == -1) {
        fatal_error(fmt::format(
          "Surface {} does not have a parameter named {} for geometry "
          "derivative {}.",
          surface_id, param, id_));
      }
      if (surface_parameter_map_.count({surface_id, param_idx}) != 0) {
        fatal_error(fmt::format("Duplicate parameter {} for surface {} in "
                                "geometry derivative {}.",
          param, surface_id, id_));
      }
      surface_parameter_map_[{surface_id, param_idx}] = n_params++;
    }
  }

  auto tally_it = model::tally_map.find(tally_id_);
  if (tally_it == model::tally_map.end()) {
    fatal_error(fmt::format(
      "Could not find tally {} for geometry derivative {}.", tally_id_, id_));
  }
  const auto& tally = model::tallies[tally_it->second];

  // ensure that the tally is using a collision estimator (for now)
  if (tally->estimator_ != TallyEstimator::COLLISION) {
    fatal_error(fmt::format(
      "Tally {} is not using a collision estimator. Geometry derivatives "
      "are only supported for collision estimators.",
      tally_id_));
  }

  // register differentiable surface parameters with surface objects
  for (const auto& [key, _] : surface_parameter_map_) {
    const auto [surface_id, param_idx] = key;
    const auto& surface = model::surfaces[model::surface_map.at(surface_id)];
    surface->register_differentiable_parameter(param_idx);
  }

  model::geometry_derivative_tallies.insert(tally_id_);

  int n_tally_bins = tally->results().shape(0) * tally->results().shape(1);
  tally_derivatives_ = tensor::Tensor<double>(
    {static_cast<size_t>(n_params), static_cast<size_t>(n_tally_bins), 10});
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
      tally_derivatives_(i, j, 7) += tally_derivatives_(i, j, 6) * norm;
      tally_derivatives_(i, j, 6) = 0.0;
      tally_derivatives_(i, j, 9) += tally_derivatives_(i, j, 8) * norm;
      tally_derivatives_(i, j, 8) = 0.0;
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
  model::geometry_derivative_tallies.clear();
  for (const auto& geom_deriv : model::geometry_derivatives) {
    geom_deriv->init_results();
    model::geometry_derivative_tallies.insert(geom_deriv->tally_id());
  }
}

void tally_geometry_derivative(int32_t geometry_deriv_idx, int64_t score_bin,
  double score, const GeometryDerivativeState& gs, const Particle& p)
{
  // std::cout << "Geom Deriv Score: " << score << std::endl;
  const auto& geom_deriv = model::geometry_derivatives[geometry_deriv_idx];
  // update the tally derivative with respect to the surface parameter
  double deriv = (gs.dj + gs.df) * score;

  int32_t data_index =
    geom_deriv->surface_parameter_map().at({gs.surface_id, gs.parameter_idx});
  // std::cout
  // << "TALLY"
  // << "\t history: " << p.id()
  // << "\t collision: " << p.n_collision()
  // << "\t cell: " << p.lowest_coord().cell()
  // << "\t material: " << p.material()
  // << "\t sigma_t: " << p.macro_xs().total
  // << "\t score: " << score
  // << "\t dj: " << gs.dj
  // << "\t df: " << gs.df
  // << "\t Sdj: " << score * gs.dj
  // << "\t Sdf: " << score * gs.df
  // << "\t total: " << score * (gs.dj + gs.df)
  // << '\n';
#pragma omp atomic
  geom_deriv->tally_derivatives()(data_index, score_bin, 0) += deriv;
#pragma omp atomic
  geom_deriv->tally_derivatives()(data_index, score_bin, 2) += score * gs.df;
#pragma omp atomic
  geom_deriv->tally_derivatives()(data_index, score_bin, 4) += score * gs.dj;
#pragma omp atomic
  geom_deriv->tally_derivatives()(data_index, score_bin, 6) += gs.df;
#pragma omp atomic
  geom_deriv->tally_derivatives()(data_index, score_bin, 8) += gs.dj;
}

void tally_geometry_derivatives(
  Particle& p, int32_t tally_id, int64_t score_bin, double score)
{
  for (const auto& geom_deriv_state : p.geometry_deriv_state()) {
    const auto& geom_deriv =
      model::geometry_derivatives[geom_deriv_state.geometry_derivative_idx];
    if (geom_deriv->tally_id() != tally_id) {
      continue;
    }
    tally_geometry_derivative(geom_deriv_state.geometry_derivative_idx,
      score_bin, score, geom_deriv_state, p);
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

extern "C" int openmc_geometry_derivative_get_surface_ids(
  int32_t index, const int32_t** ids, size_t* n)
{
  auto* deriv = get_geometry_derivative(index);
  if (!deriv)
    return OPENMC_E_OUT_OF_BOUNDS;

  *ids = deriv->surface_ids().data();
  *n = deriv->surface_ids().size();
  return 0;
}

extern "C" int openmc_geometry_derivative_get_surface_parameters(
  int32_t index, int32_t surface_id, const char** parameters)
{
  auto* deriv = get_geometry_derivative(index);
  if (!deriv)
    return OPENMC_E_OUT_OF_BOUNDS;

  for (size_t i = 0; i < deriv->surface_ids().size(); ++i) {
    if (deriv->surface_ids()[i] == surface_id) {
      *parameters = deriv->surface_parameters()[i].c_str();
      return 0;
    }
  }
  set_errmsg(fmt::format("Geometry derivative {} does not include surface {}.",
    deriv->id(), surface_id));
  return OPENMC_E_INVALID_ID;
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
