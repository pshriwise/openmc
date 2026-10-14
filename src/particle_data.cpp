#include "openmc/particle_data.h"

#include <sstream>

#include "openmc/cell.h"
#include "openmc/container_util.h"
#include "openmc/error.h"
#include "openmc/geometry.h"
#include "openmc/geometry_derivatives.h"
#include "openmc/material.h"
#include "openmc/nuclide.h"
#include "openmc/photon.h"
#include "openmc/settings.h"
#include "openmc/tallies/derivative.h"
#include "openmc/tallies/filter.h"
#include "openmc/tallies/tally.h"

namespace openmc {

void GeometryState::mark_as_lost(const char* message)
{
  fatal_error(message);
}

void GeometryState::mark_as_lost(const std::string& message)
{
  mark_as_lost(message.c_str());
}

void GeometryState::mark_as_lost(const std::stringstream& message)
{
  mark_as_lost(message.str());
}

void LocalCoord::rotate(const vector<double>& rotation)
{
  r_ = r_.rotate(rotation);
  u_ = u_.rotate(rotation);
  rotated_ = true;
}

void LocalCoord::reset()
{
  cell_ = C_NONE;
  universe_ = C_NONE;
  lattice() = C_NONE;
  lattice_index()[0] = 0;
  lattice_index()[1] = 0;
  lattice_index()[2] = 0;
  rotated() = false;
}

GeometryState::GeometryState()
{
  // Create and clear coordinate levels
  coord_.resize(model::n_coord_levels);
  cell_last_.resize(model::n_coord_levels);
  clear();
}

void GeometryState::advance_to_boundary_from_void()
{
  auto root_coord = this->coord(0);
  const auto& root_universe = model::universes[model::root_universe];
  boundary().reset();

  for (auto c_i : root_universe->cells_) {
    auto dist =
      model::cells.at(c_i)->distance(root_coord.r(), root_coord.u(), 0, this);
    if (dist.first < boundary().distance()) {
      boundary().distance() = dist.first;
      boundary().surface() = dist.second;
    }
  }

  // if no intersection or near-infinite intersection, reset
  // boundary information
  if (boundary().distance() > 1e300) {
    boundary().distance() = INFTY;
    boundary().surface() = SURFACE_NONE;
    return;
  }

  // move the particle up to (and just past) the boundary
  move_distance(boundary().distance() + TINY_BIT);
}

void GeometryState::move_distance(double length)
{
  for (int j = 0; j < n_coord(); ++j) {
    coord(j).r() += length * coord(j).u();
  }
}

ParticleData::ParticleData()
{
  zero_delayed_bank();

  // Every particle starts with no accumulated flux derivative.  Note that in
  // event mode, we construct the particle once up front, so have to run this
  // even if the current batch is inactive.
  if (!model::active_tallies.empty() || settings::event_based) {
    flux_derivs_.resize(model::tally_derivs.size());
    zero_flux_derivs();
  }

  // Allocate space for geometry derivative states
  for (int i = 0; i < model::geometry_derivatives.size(); ++i) {
    const auto& geom_deriv = model::geometry_derivatives[i];
    for (const auto& surface_parameter_pair : geom_deriv->surface_parameter_map()) {
      int surface_parameter_index = surface_parameter_pair.first;
      geometry_deriv_state().push_back(GeometryDerivativeState());
      auto& state = geometry_deriv_state().back();
      state.geometry_derivative_idx = i;
      state.surface_id = geom_deriv->surface_id();
      state.parameter_idx = surface_parameter_index;
    }
  }

  // Allocate space for tally filter matches
  filter_matches_.resize(model::tally_filters.size());

  // Create microscopic cross section caches
  neutron_xs_.resize(data::nuclides.size());
  photon_xs_.resize(data::elements.size());

  // Creates the pulse-height storage for the particle
  if (!model::pulse_height_cells.empty()) {
    pht_storage_.resize(model::pulse_height_cells.size(), 0.0);
  }
}

void ParticleData::update_geometry_derivatives()
{
  for (auto& d : geometry_deriv_state()) {
    d.update_params(*this);
  }
}


TrackState ParticleData::get_track_state() const
{
  TrackState state;
  state.r = this->r();
  state.u = this->u();
  state.E = this->E();
  state.time = this->time();
  state.wgt = this->wgt();
  state.cell_id = model::cells[this->lowest_coord().cell()]->id_;
  state.cell_instance = this->cell_instance();
  if (this->material() != MATERIAL_VOID) {
    state.material_id = model::materials[material()]->id_;
  }
  return state;
}

void GeometryDerivativeState::update_params(ParticleData& pd)
{
  if (geometry_derivative_idx == C_NONE) {
    return;
  }

  const auto& geom_deriv = model::geometry_derivatives[geometry_derivative_idx];

  // if (geom_deriv->surface_indices().count(surface_id) == 0) {
  //   return;
  // }

  const auto& surface = model::surfaces[model::surface_map[surface_id]];

  double rdot_direct = 0.0;
  // only update if we're moving toward the surface of interest
  if (surface_id == std::abs(pd.boundary().surface())) {
    rdot_direct = surface->derivative(pd.r(), pd.u(), parameter_idx);
  }

  double r = pd.boundary().distance();
  double t = std::min(pd.boundary().distance(), pd.collision_distance());

  double adv_distance = std::min(pd.boundary().distance(), pd.collision_distance());

  Position x_dot_prev = this->x_dot;

  const auto& hit_surface = model::surfaces[pd.boundary().surface_index()];

  Direction n = hit_surface->normal(pd.r() + r * pd.u());

  double rdot_inherit = -n.dot(x_dot_prev) / n.dot(pd.u());

  double rdot_total = rdot_direct + rdot_inherit;

  double tdot = rdot_total * t / r;

  double df_old = df;
  double dj_old = dj;

  f *= std::exp(-pd.macro_xs().total * adv_distance);
  df += -pd.macro_xs().total * tdot;

  this->x_dot += pd.u() * tdot;

  if (pd.boundary().distance() >= pd.collision_distance()) {
    j *= pd.boundary().distance();
    dj += rdot_total / r;
  }

  bool boundary_first = pd.boundary().distance() < pd.collision_distance();

// std::cout
//   << "GEOM"
//   << "\t history: " << pd.id()
//   << "\t collision: " << pd.n_collision()
//   << "\t event: " << (boundary_first ? "BOUNDARY" : "COLLISION")
//   << "\t r_b: " << pd.boundary().distance()
//   << "\t r_c: " << pd.collision_distance()
//   << "\t sigma_t: " << pd.macro_xs().total
//   << "\t dr_partial: " << rdot_direct
//   << "\t delta_df: " << df - df_old
//   << "\t delta_dj: " << dj - dj_old
//   << "\t df_before: " << df_old
//   << "\t df_after: " << df
//   << "\t dj_before: " << dj_old
//   << "\t dj_after: " << dj
//   << '\n';


}


} // namespace openmc
