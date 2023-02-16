#include "openmc/weight_windows.h"

#include <set>

#include "xtensor/xindex_view.hpp"
#include "xtensor/xio.hpp"
#include "xtensor/xmasked_view.hpp"
#include "xtensor/xnoalias.hpp"
#include "xtensor/xstrided_view.hpp"
#include "xtensor/xview.hpp"

#include "openmc/error.h"
#include "openmc/file_utils.h"
#include "openmc/hdf5_interface.h"
#include "openmc/mesh.h"
#include "openmc/nuclide.h"
#include "openmc/output.h"
#include "openmc/particle.h"
#include "openmc/particle_data.h"
#include "openmc/physics_common.h"
#include "openmc/search.h"
#include "openmc/tallies/filter_energy.h"
#include "openmc/tallies/filter_mesh.h"
#include "openmc/tallies/filter_particle.h"
#include "openmc/tallies/tally.h"
#include "openmc/xml_interface.h"

#include <fmt/core.h>
#include <gsl/gsl-lite.hpp>

namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

namespace variance_reduction {

std::unordered_map<int32_t, int32_t> ww_map;
openmc::vector<unique_ptr<WeightWindows>> weight_windows;

} // namespace variance_reduction

//==============================================================================
// Non-member functions
//==============================================================================

void apply_weight_windows(Particle& p)
{
  // skip dead or no energy
  if (p.E() <= 0 || !p.alive())
    return;

  bool in_domain = false;
  // TODO: this is a linear search - should do something more clever
  WeightWindow weight_window;
  for (const auto& ww : variance_reduction::weight_windows) {
    weight_window = ww->get_weight_window(p);
    if (weight_window.is_valid())
      break;
  }
  // particle is not in any of the ww domains, do nothing
  if (!weight_window.is_valid())
    return;

  // get the paramters
  double weight = p.wgt();

  // first check to see if particle should be killed for weight cutoff
  if (p.wgt() < weight_window.weight_cutoff) {
    p.wgt() = 0.0;
    return;
  }

  // check if particle is far above current weight window
  // only do this if the factor is not already set on the particle and a
  // maximum lower bound ratio is specified
  if (p.ww_factor() == 0.0 && weight_window.max_lb_ratio > 1.0 &&
      p.wgt() > weight_window.lower_weight * weight_window.max_lb_ratio) {
    p.ww_factor() =
      p.wgt() / (weight_window.lower_weight * weight_window.max_lb_ratio);
  }

  // move weight window closer to the particle weight if needed
  if (p.ww_factor() > 1.0)
    weight_window.scale(p.ww_factor());

  // if particle's weight is above the weight window split until they are within
  // the window
  if (weight > weight_window.upper_weight) {
    // do not further split the particle if above the limit
    if (p.n_split() >= settings::max_splits)
      return;

    double n_split = std::ceil(weight / weight_window.upper_weight);
    double max_split = weight_window.max_split;
    n_split = std::min(n_split, max_split);

    p.n_split() += n_split;

    // Create secondaries and divide weight among all particles
    int i_split = std::round(n_split);
    for (int l = 0; l < i_split - 1; l++) {
      p.create_secondary(weight / n_split, p.u(), p.E(), p.type());
    }
    // remaining weight is applied to current particle
    p.wgt() = weight / n_split;

  } else if (weight <= weight_window.lower_weight) {
    // if the particle weight is below the window, play Russian roulette
    double weight_survive =
      std::min(weight * weight_window.max_split, weight_window.survival_weight);
    russian_roulette(p, weight_survive);
  } // else particle is in the window, continue as normal
}

void free_memory_weight_windows()
{
  variance_reduction::ww_map.clear();
  variance_reduction::weight_windows.clear();
}

//==============================================================================
// WeightWindowSettings implementation
//==============================================================================

WeightWindows::WeightWindows(int32_t id)
{
  index_ = variance_reduction::weight_windows.size();
  set_id(id);
  set_defaults();
}

WeightWindows::WeightWindows(pugi::xml_node node)
{
  // Make sure required elements are present
  const vector<std::string> required_elems {"id", "particle_type",
    "energy_bounds", "lower_ww_bounds", "upper_ww_bounds"};
  for (const auto& elem : required_elems) {
    if (!check_for_node(node, elem.c_str())) {
      fatal_error(fmt::format("Must specify <{}> for weight windows.", elem));
    }
  }

  // Get weight windows ID
  int32_t id = std::stoi(get_node_value(node, "id"));
  this->set_id(id);

  // get the particle type
  auto particle_type_str = std::string(get_node_value(node, "particle_type"));
  particle_type_ = openmc::str_to_particle_type(particle_type_str);

  // Determine associated mesh
  int32_t mesh_id = std::stoi(get_node_value(node, "mesh"));
  mesh_idx_ = model::mesh_map.at(mesh_id);

  // energy bounds
  if (check_for_node(node, "energy_bounds"))
    energy_bounds_ = get_node_array<double>(node, "energy_bounds");

  // get the survival value - optional
  if (check_for_node(node, "survival_ratio")) {
    survival_ratio_ = std::stod(get_node_value(node, "survival_ratio"));
    if (survival_ratio_ <= 1)
      fatal_error("Survival to lower weight window ratio must bigger than 1 "
                  "and less than the upper to lower weight window ratio.");
  }

  // get the max lower bound ratio - optional
  if (check_for_node(node, "max_lower_bound_ratio")) {
    max_lb_ratio_ = std::stod(get_node_value(node, "max_lower_bound_ratio"));
    if (max_lb_ratio_ < 1.0) {
      fatal_error("Maximum lower bound ratio must be larger than 1");
    }
  }

  // get the max split - optional
  if (check_for_node(node, "max_split")) {
    max_split_ = std::stod(get_node_value(node, "max_split"));
    if (max_split_ <= 1)
      fatal_error("max split must be larger than 1");
  }

  // weight cutoff - optional
  if (check_for_node(node, "weight_cutoff")) {
    weight_cutoff_ = std::stod(get_node_value(node, "weight_cutoff"));
    if (weight_cutoff_ <= 0)
      fatal_error("weight_cutoff must be larger than 0");
    if (weight_cutoff_ > 1)
      fatal_error("weight_cutoff must be less than 1");
  }

  // read the lower/upper weight bounds
  this->set_weight_windows(get_node_array<double>(node, "lower_ww_bounds"),
    get_node_array<double>(node, "upper_ww_bounds"));

  set_defaults();
}

WeightWindows::~WeightWindows()
{
  variance_reduction::ww_map.erase(id());
}

WeightWindows* WeightWindows::create(int32_t id)
{
  variance_reduction::weight_windows.push_back(make_unique<WeightWindows>());
  auto wws = variance_reduction::weight_windows.back().get();
  variance_reduction::ww_map[wws->id()] =
    variance_reduction::weight_windows.size() - 1;
  return wws;
}

WeightWindows* WeightWindows::from_hdf5(
  hid_t wws_group, const std::string& group_name)
{
  // collect ID from the name of this group
  hid_t ww_group = open_group(wws_group, group_name.c_str());

  auto wws = WeightWindows::create();

  std::string particle_type;
  read_dataset(ww_group, "particle_type", particle_type);
  wws->particle_type_ = openmc::str_to_particle_type(particle_type);

  read_dataset<double>(ww_group, "energy_bounds", wws->energy_bounds_);

  int32_t mesh_id;
  read_dataset(ww_group, "mesh", mesh_id);

  if (model::mesh_map.count(mesh_id) == 0) {
    fatal_error(
      fmt::format("Mesh {} used in weight windows does not exist.", mesh_id));
  }
  wws->set_mesh(model::mesh_map[mesh_id]);

  wws->lower_ww_ = xt::empty<double>(wws->bounds_size());
  wws->upper_ww_ = xt::empty<double>(wws->bounds_size());

  read_dataset<double>(ww_group, "lower_ww_bounds", wws->lower_ww_);
  read_dataset<double>(ww_group, "upper_ww_bounds", wws->upper_ww_);
  read_dataset(ww_group, "survival_ratio", wws->survival_ratio_);
  read_dataset(ww_group, "max_lower_bound_ratio", wws->max_lb_ratio_);
  read_dataset(ww_group, "max_split", wws->max_split_);
  read_dataset(ww_group, "weight_cutoff", wws->weight_cutoff_);

  close_group(ww_group);

  return wws;
}

void WeightWindows::set_defaults()
{
  // ensure default values are set
  if (energy_bounds_.size() == 0) {
    int p_type = static_cast<int>(particle_type_);
    energy_bounds_.push_back(data::energy_min[p_type]);
    energy_bounds_.push_back(data::energy_max[p_type]);
  }
}

void WeightWindows::set_id(int32_t id)
{
  Expects(id >= 0 || id == C_NONE);

  // Clear entry in mesh map in case one was already assigned
  if (id_ != C_NONE) {
    variance_reduction::ww_map.erase(id_);
    id_ = C_NONE;
  }

  // Ensure no other mesh has the same ID
  if (variance_reduction::ww_map.find(id) != variance_reduction::ww_map.end()) {
    throw std::runtime_error {
      fmt::format("Two weight windows have the same ID: {}", id)};
  }

  // If no ID is specified, auto-assign the next ID in the sequence
  if (id == C_NONE) {
    id = 0;
    for (const auto& m : variance_reduction::weight_windows) {
      id = std::max(id, m->id_);
    }
    ++id;
  }

  // Update ID and entry in the mesh map
  id_ = id;
  variance_reduction::ww_map[id] = index_;
}

void WeightWindows::set_energy_bounds(gsl::span<double const> bounds)
{
  energy_bounds_.clear();
  energy_bounds_.insert(energy_bounds_.begin(), bounds.begin(), bounds.end());
  // TODO: check that sizes still make sense
}

void WeightWindows::set_particle_type(ParticleType p_type)
{
  if (p_type != ParticleType::neutron && p_type != ParticleType::photon)
    fatal_error(
      fmt::format("Particle type '{}' cannot be applied to weight windows.",
        particle_type_to_str(p_type)));
  particle_type_ = p_type;
}

void WeightWindows::set_mesh(int32_t mesh_idx)
{
  if (mesh_idx < 0 || mesh_idx >= model::meshes.size())
    fatal_error(fmt::format("Could not find a mesh for index {}", mesh_idx));

  mesh_idx_ = mesh_idx;
}

void WeightWindows::set_mesh(const std::unique_ptr<Mesh>& mesh)
{
  set_mesh(mesh.get());
}

void WeightWindows::set_mesh(const Mesh* mesh)
{
  set_mesh(model::mesh_map[mesh->id_]);
}

WeightWindow WeightWindows::get_weight_window(const Particle& p) const
{
  // check for particle type
  if (particle_type_ != p.type()) {
    return {};
  }

  // Get mesh index for particle's position
  const auto& mesh = this->mesh();
  int mesh_bin = mesh->get_bin(p.r());

  // particle is outside the weight window mesh
  if (mesh_bin < 0)
    return {};

  // particle energy
  double E = p.E();

  // check to make sure energy is in range, expects sorted energy values
  if (E < energy_bounds_.front() || E > energy_bounds_.back())
    return {};

  // get the mesh bin in energy group
  int energy_bin =
    lower_bound_index(energy_bounds_.begin(), energy_bounds_.end(), E);

  // mesh_bin += energy_bin * mesh->n_bins();
  // Create individual weight window
  WeightWindow ww;
  ww.lower_weight = lower_ww_(energy_bin, mesh_bin);
  ww.upper_weight = upper_ww_(energy_bin, mesh_bin);
  ww.survival_weight = ww.lower_weight * survival_ratio_;
  ww.max_lb_ratio = max_lb_ratio_;
  ww.max_split = max_split_;
  ww.weight_cutoff = weight_cutoff_;
  return ww;
}

std::array<int, 2> WeightWindows::bounds_size() const
{
  int num_spatial_bins = this->mesh()->n_bins();
  int num_energy_bins =
    energy_bounds_.size() > 0 ? energy_bounds_.size() - 1 : 1;
  return {num_energy_bins, num_spatial_bins};
}

template<class T>
void WeightWindows::check_bounds(const T& lower, const T& upper) const
{
  // make sure that the upper and lower bounds have the same size
  if (lower.size() != upper.size()) {
    auto msg = fmt::format("The upper and lower weight window lengths do not "
                           "match.\n Lower size: {}\n Upper size: {}",
      lower.size(), upper.size());
    fatal_error(msg);
  }
  check_bounds(lower);
}

template<class T>
void WeightWindows::check_bounds(const T& bounds) const
{
  // check that the number of weight window entries is correct
  auto dims = this->bounds_size();
  if (bounds.size() != dims[0] * dims[1]) {
    auto err_msg =
      fmt::format("In weight window domain {} the number of spatial "
                  "energy/spatial bins ({}) does not match the number "
                  "of weight bins ({})",
        id_, this->bounds_size(), bounds.size());
    fatal_error(err_msg);
  }
}

void WeightWindows::set_weight_windows(
  const xt::xtensor<double, 2>& lower_bounds,
  const xt::xtensor<double, 2>& upper_bounds)
{

  check_bounds(lower_bounds, upper_bounds);

  // set new weight window values
  lower_ww_ = lower_bounds;
  upper_ww_ = upper_bounds;
}

void WeightWindows::set_weight_windows(
  const xt::xtensor<double, 2>& lower_bounds, double ratio)
{
  check_bounds(lower_bounds);

  // set new weight window values
  lower_ww_ = lower_bounds;
  upper_ww_ = lower_bounds;
  upper_ww_ *= ratio;
}

void WeightWindows::set_weight_windows(
  gsl::span<const double> lower_bounds, gsl::span<const double> upper_bounds)
{
  check_bounds(lower_bounds, upper_bounds);
  lower_ww_ = xt::empty<double>(bounds_size());
  upper_ww_ = xt::empty<double>(bounds_size());

  // set new weight window values
  auto shape = bounds_size();
  xt::view(lower_ww_, xt::all()) = xt::adapt(lower_bounds.data(), lower_ww_.shape());
  xt::view(upper_ww_, xt::all()) = xt::adapt(upper_bounds.data(), upper_ww_.shape());
}

void WeightWindows::set_weight_windows(
  gsl::span<const double> lower_bounds, double ratio)
{
  check_bounds(lower_bounds);

  lower_ww_ = xt::empty<double>(bounds_size());
  upper_ww_ = xt::empty<double>(bounds_size());

  // set new weight window values
  auto shape = bounds_size();
  xt::view(lower_ww_, xt::all()) = xt::adapt(lower_bounds.data(), lower_ww_.shape());
  xt::view(upper_ww_, xt::all()) = xt::adapt(lower_bounds.data(), upper_ww_.shape());
  upper_ww_ *= ratio;
}

void WeightWindows::update_weight_windows_magic(
  const Tally* tally, const std::string& value, double threshold, double ratio)
{

  ///////////////////////////
  // Setup and checks
  ///////////////////////////
  check_tally_update_compatibility(tally);

  auto filter_indices = tally->filter_indices();

  // empty out previous results for weight window bounds
  if (lower_ww_.size() == 0 || upper_ww_.size() == 0) {
    lower_ww_ = xt::empty<double>(bounds_size());
    upper_ww_ = xt::empty<double>(bounds_size());
  }

  // adjust filter indices for dummy axes we may introduce when reshaping tally data
  if (!filter_indices.count(FilterType::PARTICLE)) {
    if (!filter_indices.count(FilterType::ENERGY)) {
      filter_indices[FilterType::MESH] += 2;
    } else {
      filter_indices[FilterType::MESH] += 1;
      filter_indices[FilterType::ENERGY] += 1;
    }
  }

  // sanitize filter indices
  for (auto& pair : filter_indices) {
    if (pair.second > 3) pair.second -= 3;
  }

  // determine which value to use
  const std::set<std::string> allowed_values = {"mean", "rel_err"};
  if (allowed_values.count(value) == 0) {
    fatal_error(fmt::format("Invalid value '{}' specified for weight window "
                            "generation. Must be one of: 'mean' or 'rel_err'",
      value));
  }

  // determine the index of the specified score
  int score_index = tally->score_index("flux");
  if (score_index == C_NONE) {
    fatal_error(
      fmt::format("A 'flux' score required for weight window generation "
                  "is not present on tally {}.",
        tally->id()));
  }

  ///////////////////////////
  // Extract tally data
  //
  // At the end of this section, the mean and rel_err array
  // is a 2D view of tally data (n_mesh_bins, n_e_groups)
  //
  ///////////////////////////

  // build a shape for a view of the tally results, this will always be
  // dimension 5 (3 filter dimensions, 1 score dimension, 1 results dimension)
  std::array<int, 5> shape = {1, 1, 1, tally->n_scores(), static_cast<int>(TallyResult::TALLY_RESULT_SIZE)};
  // build the transpose information to re-order data according to filter type
  std::array<int, 5> transpose = {0, 1, 2, 3, 4};

  // determine the dimension and index of the particle data
  int particle_idx = 0;
  if (filter_indices.count(FilterType::PARTICLE)) {
    // get the particle filter
    auto pf = tally->get_filter<ParticleFilter>();
    const auto& particles = pf->particles();

    // find the index of the particle that matches these weight windows
    auto p_it =
      std::find(particles.begin(), particles.end(), this->particle_type_);
    // if the particle filter doesn't have particle data for the particle
    // used on this weight windows instance, report an error
    if (p_it == particles.end()) {
      auto msg = fmt::format("Particle type '{}' not present on Filter {} for "
                             "Tally {} used to update WeightWindows {}",
        particle_type_to_str(this->particle_type_), pf->id(), tally->id(),
        this->id());
      fatal_error(msg);
    }

    // use the index of the particle in the filter to down-select data later
    particle_idx = p_it - particles.begin();
    shape[filter_indices[FilterType::PARTICLE]] = pf->n_bins();
    transpose[0] = filter_indices[FilterType::PARTICLE];
  }

  if (filter_indices.count(FilterType::ENERGY)) {
    auto ef = tally->get_filter<EnergyFilter>();
    shape[filter_indices[FilterType::ENERGY]] = ef->n_bins();
    transpose[1] = filter_indices[FilterType::ENERGY];
  }

  if (filter_indices.count(FilterType::MESH)) {
    auto mf = tally->get_filter<MeshFilter>();
    shape[filter_indices[FilterType::MESH]] = mf->n_bins();
    transpose[2] = filter_indices[FilterType::MESH];
  }

  // get a fully reshaped view of the tally according to tally ordering of filters
  auto tally_values = xt::reshape_view(tally->results(), shape);

  // get a that is (particle, energy, mesh, scores, values)
  auto transposed_view = xt::transpose(tally_values, transpose);

  // down-select data based on particle and score
  auto sum = xt::view(transposed_view, particle_idx, xt::all(), xt::all(),
    score_index, static_cast<int>(TallyResult::SUM));
  auto sum_sq = xt::view(transposed_view, particle_idx, xt::all(), xt::all(), score_index, static_cast<int>(TallyResult::SUM_SQ));
  int n = tally->n_realizations_;

  //////////////////////////////////////////////
  //
  // Assign new weight windows
  //
  // Use references to the existing weight window data
  // to store and update the values
  //
  //////////////////////////////////////////////

  // up to this point the data arrays are views into the tally results (no
  // computation has been performed) now we'll switch references to the tally's
  // bounds to avoid allocating additional memory
  auto& new_bounds  = this->lower_ww_;
  // noalias avoids memory allocation here
  xt::noalias(new_bounds) = sum / n;

  auto& rel_err = this->upper_ww_;
  rel_err =
    xt::sqrt(((sum_sq / n) - xt::square(new_bounds)) / (n - 1)) / new_bounds;
  xt::filter(rel_err, sum <= 0.0).fill(INFTY);

  if (value == "rel_err")
    new_bounds = 1 / rel_err;

  // get mesh volumes
  auto mesh_vols = this->mesh()->volumes();

  int e_bins = new_bounds.shape()[0];
  for (int e = 0; e < e_bins; e++) {
    // select all
    auto group_view = xt::view(new_bounds, e);

    double group_max = *std::max_element(group_view.begin(), group_view.end());
    // normalize values in this energy group by the maximum value for this
    // group
    if (group_max > 0.0)
      group_view /= group_max;

    // divide by volume of mesh elements
    for (int i = 0; i < group_view.size(); i++) {
      group_view[i] /= mesh_vols[i];
    }
  }

  // make sure that values where the mean is zero are set s.t. the weight window
  // value will be ignored
  xt::filter(new_bounds, sum <= 0.0).fill(-1.0);

  // make sure the weight windows are ignored for any locations where the
  // relative error is higher than the specified relative error threshold
  xt::filter(new_bounds, rel_err > threshold).fill(-1.0);

  // update the bounds of this weight window class
  // noalias avoids additional memory allocation
  xt::noalias(upper_ww_) = ratio * lower_ww_;
}

void WeightWindows::check_tally_update_compatibility(const Tally* tally)
{
  // define the set of allowed filters for the tally
  const std::set<FilterType> allowed_filters = {
    FilterType::MESH, FilterType::ENERGY, FilterType::PARTICLE};

  // retrieve a mapping of filter type to filter index for the tally
  auto filter_indices = tally->filter_indices();

  // a mesh filter is required for a tally used to update weight windows
  if (!filter_indices.count(FilterType::MESH)) {
    fatal_error(
      "A mesh filter is required for a tally to update weight window bounds");
  }

  // ensure the mesh filter is using the same mesh as this weight window object
  auto mesh_filter = tally->get_filter<MeshFilter>();

  // make sure that all of the filters present on the tally are allowed
  for (auto filter_pair : filter_indices) {
    if (allowed_filters.find(filter_pair.first) == allowed_filters.end()) {
      fatal_error(fmt::format("Invalid filter type '{}' found on tally "
                              "used for weight window generation.",
        model::tally_filters[filter_pair.second]->type_str()));
    }
  }

  if (mesh_filter->mesh() != mesh_idx_) {
    int32_t mesh_filter_id = model::meshes[mesh_filter->mesh()]->id();
    int32_t ww_mesh_id = model::meshes[this->mesh_idx_]->id();
    fatal_error(fmt::format("Mesh filter {} uses a different mesh ({}) than "
                            "weight window {} mesh ({})",
      mesh_filter->id(), mesh_filter_id, id_, ww_mesh_id));
  }

  // if an energy filter exists, make sure the energy grid matches that of this
  // weight window object
  if (auto energy_filter = tally->get_filter<EnergyFilter>()) {
    std::vector<double> filter_bins = energy_filter->bins();
    std::set<double> filter_e_bounds(
      energy_filter->bins().begin(), energy_filter->bins().end());
    if (filter_e_bounds.size() != energy_bounds().size()) {
      fatal_error(
        fmt::format("Energy filter {} does not have the same number of energy "
                    "bounds ({}) as weight window object {} ({})",
          energy_filter->id(), filter_e_bounds.size(), id_,
          energy_bounds().size()));
    }

    for (auto e : energy_bounds()) {
      if (filter_e_bounds.count(e) == 0) {
        fatal_error(fmt::format(
          "Energy bounds of filter {} and weight windows {} do not match",
          energy_filter->id(), id_));
      }
    }
  }
}

void WeightWindows::export_to_hdf5(const std::string& filename) const
{
  hid_t file_id = file_open(filename, 'w');

  // Write file type
  write_attribute(file_id, "filetype", "weight_windows");

  // Write revisiion number for state point file
  write_attribute(file_id, "version", VERSION_WEIGHT_WINDOWS);

  hid_t weight_windows_group =
    create_group(file_id, fmt::format("weight_windows", id_));

  // write instance information to file
  this->to_hdf5(weight_windows_group);

  close_group(weight_windows_group);

  file_close(file_id);
}

void WeightWindows::to_hdf5(hid_t group) const
{
  hid_t ww_group = create_group(group, fmt::format("weight_windows {}", id()));

  write_dataset(ww_group, "mesh", this->mesh()->id());
  write_dataset(
    ww_group, "particle_type", openmc::particle_type_to_str(particle_type_));
  write_dataset(ww_group, "energy_bounds", energy_bounds_);
  write_dataset(ww_group, "lower_ww_bounds", lower_ww_);
  write_dataset(ww_group, "upper_ww_bounds", upper_ww_);
  write_dataset(ww_group, "survival_ratio", survival_ratio_);
  write_dataset(ww_group, "max_lower_bound_ratio", max_lb_ratio_);
  write_dataset(ww_group, "max_split", max_split_);
  write_dataset(ww_group, "weight_cutoff", weight_cutoff_);

  close_group(ww_group);
}

//==============================================================================
// C API
//==============================================================================

int verify_ww_index(int32_t index)
{
  if (index < 0 || index >= variance_reduction::weight_windows.size()) {
    set_errmsg(fmt::format("Index '{}' for weight windows is invalid", index));
    return OPENMC_E_OUT_OF_BOUNDS;
  }
  return 0;
}

extern "C" int openmc_weight_windows_get_index(int32_t id, int32_t* idx)
{
  auto it = variance_reduction::ww_map.find(id);
  if (it == variance_reduction::ww_map.end()) {
    set_errmsg(fmt::format("No weight windows exist with ID={}", id));
    return OPENMC_E_INVALID_ID;
  }

  *idx = it->second;
  return 0;
}

extern "C" int openmc_weight_windows_get_id(int32_t index, int32_t* id)
{
  if (int err = verify_ww_index(index))
    return err;

  const auto& wws = variance_reduction::weight_windows.at(index);
  *id = wws->id();
  return 0;
}

extern "C" int openmc_weight_windows_set_id(int32_t index, int32_t id)
{
  if (int err = verify_ww_index(index))
    return err;

  const auto& wws = variance_reduction::weight_windows.at(index);
  wws->id() = id;
  variance_reduction::ww_map[id] = index;
  return 0;
}

extern "C" int openmc_set_weight_windows(
  int ww_id, size_t n, const double* lower_bounds, const double* upper_bounds)
{

  // look up the weight windows object
  const auto& wws =
    variance_reduction::weight_windows.at(variance_reduction::ww_map.at(ww_id));

  // check length of arrays
  auto dims = wws->bounds_size();
  if (n != dims[0] * dims[1]) {
    set_errmsg(fmt::format(
      "Incorrect size for weight window bounds for domain {}", wws->id()));
    return OPENMC_E_INVALID_ARGUMENT;
  }

  // set bounds
  wws->set_weight_windows({lower_bounds, n}, {upper_bounds, n});

  return 0;
}

extern "C" int openmc_update_weight_windows_magic(int32_t tally_idx,
  int32_t ww_idx, const char* value, double threshold, double ratio)
{
  // get the requested tally
  const Tally* tally = model::tallies.at(tally_idx).get();

  // get the WeightWindows object
  const auto& wws = variance_reduction::weight_windows.at(ww_idx);

  wws->update_weight_windows_magic(tally, value, threshold, ratio);

  return 0;
}

extern "C" int openmc_weight_windows_set_mesh(int32_t ww_idx, int32_t mesh_idx)
{
  const auto& wws = variance_reduction::weight_windows.at(ww_idx);
  wws->set_mesh(mesh_idx);
  return 0;
}

extern "C" int openmc_weight_windows_get_mesh(int32_t ww_idx, int32_t* mesh_idx)
{
  const auto& wws = variance_reduction::weight_windows.at(ww_idx);
  *mesh_idx = model::mesh_map.at(wws->mesh()->id());
  return 0;
}

extern "C" int openmc_weight_windows_set_energy_bounds(
  int32_t ww_idx, double* e_bounds, size_t e_bounds_size)
{
  const auto& wws = variance_reduction::weight_windows.at(ww_idx);
  wws->set_energy_bounds({e_bounds, e_bounds_size});
  return 0;
}

extern "C" int openmc_weight_windows_get_energy_bounds(
  int32_t ww_idx, const double** e_bounds, size_t* e_bounds_size)
{

  const auto& wws = variance_reduction::weight_windows[ww_idx].get();
  *e_bounds = wws->energy_bounds().data();
  *e_bounds_size = wws->energy_bounds().size();
  return 0;
}

extern "C" int openmc_weight_windows_set_particle(
  int32_t index, const char* particle)
{
  if (int err = verify_ww_index(index))
    return err;

  const auto& wws = variance_reduction::weight_windows.at(index);
  wws->set_particle_type(str_to_particle_type(particle));
  return 0;
}

extern "C" int openmc_weight_windows_get_bounds(int32_t index,
  const double** lower_bounds, const double** upper_bounds, size_t* size)
{
  if (int err = verify_ww_index(index))
    return err;

  const auto& wws = variance_reduction::weight_windows[index];
  *size = wws->lower_ww_bounds().size();
  *lower_bounds = wws->lower_ww_bounds().data();
  *upper_bounds = wws->upper_ww_bounds().data();
  return 0;
}

extern "C" int openmc_weight_windows_set_bounds(int32_t index,
  const double* lower_bounds, const double* upper_bounds, size_t size)
{
  if (int err = verify_ww_index(index))
    return err;

  const auto& wws = variance_reduction::weight_windows[index];
  wws->set_weight_windows({lower_bounds, size}, {upper_bounds, size});
  return 0;
}

extern "C" int openmc_weight_windows_get_particle(int32_t index, int* particle)
{
  if (int err = verify_ww_index(index))
    return err;

  const auto& wws = variance_reduction::weight_windows.at(index);
  *particle = static_cast<int>(wws->particle_type());
  return 0;
}

extern "C" int openmc_extend_weight_windows(
  int32_t n, int32_t* index_start, int32_t* index_end)
{
  if (index_start)
    *index_start = variance_reduction::weight_windows.size();
  if (index_end)
    *index_end = variance_reduction::weight_windows.size() + n - 1;
  for (int i = 0; i < n; ++i)
    variance_reduction::weight_windows.push_back(
      make_unique<WeightWindows>(-1));
  return 0;
}

extern "C" size_t openmc_weight_windows_size()
{
  return variance_reduction::weight_windows.size();
}

extern "C" int openmc_weight_windows_export(const char* filename)
{

  std::string name = filename ? filename : "weight_windows.h5";

  write_message(fmt::format("Exporting weight windows to {}...", name), 5);

  hid_t ww_file = file_open(name, 'w');

  // Write file type
  write_attribute(ww_file, "filetype", "weight_windows");

  // Write revisiion number for state point file
  write_attribute(ww_file, "version", VERSION_WEIGHT_WINDOWS);

  hid_t weight_windows_group = create_group(ww_file, "weight_windows");

  hid_t mesh_group = create_group(ww_file, "meshes");

  for (const auto& ww : variance_reduction::weight_windows) {
    ww->to_hdf5(weight_windows_group);
    ww->mesh()->to_hdf5(mesh_group);
  }

  close_group(weight_windows_group);

  file_close(ww_file);

  return 0;
}

extern "C" int openmc_weight_windows_import(const char* filename)
{

  std::string name = filename ? filename : "weight_windows.h5";

  write_message(fmt::format("Importing weight windows from {}...", name), 5);

  if (!file_exists(name)) {
    set_errmsg(fmt::format("File '{}' does not exist", name));
  }

  hid_t ww_file = file_open(name, 'r');

  // Check that filetype is correct
  std::string filetype;
  read_attribute(ww_file, "filetype", filetype);
  if (filetype != "weight_windows") {
    file_close(ww_file);
    set_errmsg(fmt::format("File '{}' is not a weight windows file.", name));
    return OPENMC_E_INVALID_ARGUMENT;
  }

  // Check that the file version is compatible
  std::array<int, 2> file_version;
  read_attribute(ww_file, "version", file_version);
  if (file_version[0] != VERSION_WEIGHT_WINDOWS[0]) {
    std::string err_msg =
      fmt::format("File '{}' has version {} which is incompatible with the "
                  "expected version ({}).",
        name, file_version, VERSION_WEIGHT_WINDOWS);
    set_errmsg(err_msg);
    return OPENMC_E_INVALID_ARGUMENT;
  }

  hid_t weight_windows_group = open_group(ww_file, "weight_windows");

  std::vector<std::string> names = group_names(weight_windows_group);

  for (const auto& name : names) {
    WeightWindows::from_hdf5(weight_windows_group, name);
  }

  close_group(weight_windows_group);

  file_close(ww_file);

  return 0;
}

} // namespace openmc
