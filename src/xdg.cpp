#include "openmc/xdg.h"

#include "openmc/constants.h"
#include "openmc/container_util.h"
#include "openmc/error.h"
#include "openmc/file_utils.h"
#include "openmc/geometry.h"
#include "openmc/geometry_aux.h"
#include "openmc/hdf5_interface.h"
#include "openmc/material.h"
#include "openmc/settings.h"
#include "openmc/string_utils.h"

#include <fmt/core.h>

#ifdef OPENMC_XDG
#include "xdg/xdg.h"
#endif

#include <algorithm>
#include <fstream>
#include <sstream>
#include <string>

namespace openmc {

#ifdef OPENMC_XDG
const bool XDG_ENABLED = true;
#else
const bool XDG_ENABLED = false;
#endif

} // namespace openmc

#ifdef OPENMC_XDG

namespace openmc {

//==============================================================================
// XDG Universe implementation
//==============================================================================

XDGUniverse::XDGUniverse(pugi::xml_node node)
{
  if (check_for_node(node, "id")) {
    id_ = std::stoi(get_node_value(node, "id"));
  } else {
    fatal_error("Must specify the id of the XDG universe");
  }

  if (check_for_node(node, "filename")) {
    filename_ = get_node_value(node, "filename");
    if (!starts_with(filename_, "/")) {
      filename_ = dir_name(settings::path_input) + filename_;
    }
  } else {
    fatal_error("Must specify a file for the XDG universe");
  }

  adjust_geometry_ids_ = false;
  if (check_for_node(node, "auto_geom_ids")) {
    adjust_geometry_ids_ = get_node_value_bool(node, "auto_geom_ids");
  }

  adjust_material_ids_ = false;
  if (check_for_node(node, "auto_mat_ids")) {
    adjust_material_ids_ = get_node_value_bool(node, "auto_mat_ids");
  }

  initialize();
}

XDGUniverse::XDGUniverse(
  const std::string& filename, bool auto_geom_ids, bool auto_mat_ids)
  : filename_(filename), adjust_geometry_ids_(auto_geom_ids),
    adjust_material_ids_(auto_mat_ids)
{
  geom_type() = GeometryType::XDG;
  set_id();
  initialize();
}

XDGUniverse::XDGUniverse(std::shared_ptr<xdg::XDG> xdg_ptr,
  const std::string& filename, bool auto_geom_ids, bool auto_mat_ids)
  : xdg_instance_(xdg_ptr), filename_(filename),
    adjust_geometry_ids_(auto_geom_ids), adjust_material_ids_(auto_mat_ids)
{
  geom_type() = GeometryType::XDG;
  set_id();
  init_metadata();
  init_geometry();
}

void XDGUniverse::set_id()
{
  // determine the next universe id
  int32_t next_univ_id = 0;
  for (const auto& u : model::universes) {
    if (u->id_ > next_univ_id)
      next_univ_id = u->id_;
  }
  next_univ_id++;

  // set the universe id
  id_ = next_univ_id;
}

void XDGUniverse::initialize()
{
  geom_type() = GeometryType::XDG;

  init_xdg();
  init_metadata();
  init_geometry();
}

void XDGUniverse::init_xdg()
{
  // create a new XDG instance
  xdg_instance_ = xdg::XDG::create(xdg::MeshLibrary::MOAB);

  // load the XDG geometry
  if (!file_exists(filename_)) {
    fatal_error("Geometry XDG file '" + filename_ + "' does not exist!");
  }
  xdg_instance_->mesh_manager()->load_file(filename_);

  xdg_instance_->mesh_manager()->init();
  xdg_instance_->mesh_manager()->parse_metadata();
  xdg_instance_->prepare_raytracer();
}

void XDGUniverse::init_metadata()
{
  xdg_instance_->mesh_manager()->parse_metadata();
}

void XDGUniverse::init_geometry()
{
  int32_t next_cell_id = 0;
  for (const auto& c : model::cells) {
    if (c->id_ > next_cell_id)
      next_cell_id = c->id_;
  }
  cell_idx_offset_ = model::cells.size();
  next_cell_id++;

  for (auto volume : xdg_ptr()->mesh_manager()->volumes()) {
    auto c = std::make_unique<XDGCell>(xdg_instance_, volume);
    c->id_ =     c->id_ = adjust_geometry_ids_ ? next_cell_id++ : volume;
    c->universe_ = this->id_;

    auto in_map = model::cell_map.find(c->id_);
    if (in_map == model::cell_map.end()) {
      model::cell_map[c->id_] = model::cells.size();
      cell_index_map_[volume] = model::cells.size();
    } else {
      warning(fmt::format("XDG Cell IDs: {}", xdg_ids_for_dim(3)));
      fatal_error(fmt::format(
        "XDG Universe {} contains a cell with ID {}, which "
        "already exists elsewhere in the geometry. Setting auto_geom_ids "
        "to True when defining the XDG Universe may "
        "resolve this issue",
        this->id_, c->id_));
    }

    // --- Materials ---

    // determine volume material assignment
    std::string mat_str;
    if (!xdg_ptr()->mesh_manager()->volume_has_property(volume, xdg::PropertyType::MATERIAL)) {
      fatal_error(fmt::format("XDG Volume {} has no material assignment."), c->id_);
    } else {
      mat_str = xdg_ptr()->mesh_manager()->get_volume_property(volume, xdg::PropertyType::MATERIAL).value;
    }

    to_lower(mat_str);

    assign_material(mat_str, c);

    // ----Temperatures-----

    // if the temperature property is set on the XDG volume, use that
    if (xdg_ptr()->mesh_manager()->volume_has_property(volume, xdg::PropertyType::TEMPERATURE)) {
      // read the temperature from the file (assumed to be in K)
      xdg::Property xdg_temperature = xdg_ptr()->mesh_manager()->get_volume_property(volume, xdg::PropertyType::TEMPERATURE);
      double temperature = std::stod(xdg_temperature.value);
      c->sqrtkT_.push_back(std::sqrt(K_BOLTZMANN * temperature));
    } else { // otherwise, set the temperature as usual
      const auto& mat = model::materials[model::material_map[c->material_[0]]];
      if (mat->temperature() > 0.0) {
        c->sqrtkT_.push_back(std::sqrt(K_BOLTZMANN * mat->temperature()));
      } else {
      c->sqrtkT_.push_back(
        std::sqrt(K_BOLTZMANN * settings::temperature_default));
      }
    }
    model::cells.emplace_back(std::move(c));
  }

  // allocate the cell overlap count if necessary
  if (settings::check_overlaps) {
    model::overlap_check_count.resize(model::cells.size(), 0);
  }
  // determine the next surface id
  int32_t next_surf_id = 0;
  for (const auto& s : model::surfaces) {
    if (s->id_ > next_surf_id)
      next_surf_id = s->id_;
  }
  surf_idx_offset_ = model::surfaces.size();
  next_surf_id++;

  for (auto surface : xdg_ptr()->mesh_manager()->surfaces()) {
    auto s = std::make_unique<XDGSurface>(xdg_instance_, surface);
    s->id_ = adjust_geometry_ids_ ? next_surf_id++ : surface;

    s->surf_source_ = contains(settings::source_write_surf_id, s->id_);

    // ----- Boundary Conditions--------
    if (xdg_ptr()->mesh_manager()->surface_has_property(surface, xdg::PropertyType::BOUNDARY_CONDITION)) {
      xdg::Property xdg_bc = xdg_ptr()->mesh_manager()->get_surface_property(surface, xdg::PropertyType::BOUNDARY_CONDITION);
      std::string bc_value = xdg_bc.value;
      to_lower(bc_value);

      if (bc_value.empty() || bc_value == "transmit" ||
          bc_value == "transmission") {
        // set to transmission by default (nullptr)
      } else if (bc_value == "vacuum") {
        s->bc_ = make_unique<VacuumBC>();
      } else if (bc_value == "reflective" || bc_value == "reflect" ||
                bc_value == "reflecting") {
        s->bc_ = make_unique<ReflectiveBC>();
      } else if (bc_value == "periodic") {
        fatal_error("Periodic boundary condition not supported in XDG geometry.");
      } else {
        fatal_error(fmt::format("Unknown boundary condition \"{}\" specified "
                                "on surface {}",
          bc_value, s->id_));
      }
    }

    // add to global array and map
    auto in_map = model::surface_map.find(s->id_);
    if (in_map == model::surface_map.end()) {
      model::surface_map[s->id_] = model::surfaces.size();
      surface_index_map_[surface] = model::surfaces.size() + 1;
    } else {
      warning(fmt::format("XDG Surface IDs: {}", xdg_ids_for_dim(2)));
      fatal_error(fmt::format("Surface ID {} exists in both Universe {} "
                              "and the CSG geometry.",
        s->id_, this->id_));
    }

    model::surfaces.emplace_back(std::move(s));
  }
}

int32_t XDGUniverse::cell_index(xdg::MeshID volume) const
{
  // return the index of the volume in the XDG instance and then
  // adjust by the offset into the model cells for this XDG universe
  return cell_index_map_.at(volume);
}

int32_t XDGUniverse::surface_index(xdg::MeshID surface) const
{
  // return the index of the surface in the XDG instance and then
  // adjust by the offset into the model cells for this XDG universe
  return surface_index_map_.at(surface);
}

std::string XDGUniverse::xdg_ids_for_dim(int dim) const
{
  // get the ID vector
  std::vector<xdg::MeshID> id_vec;
  if (dim == 3) {
    id_vec = xdg_ptr()->mesh_manager()->volumes();
  } else if (dim == 2) {
    id_vec = xdg_ptr()->mesh_manager()->surfaces();
  } else {
    warning(fmt::format("Invalid dimension {} for XDG ID generation", dim));
    return "";
  }

  // sort the vector of ids
  std::sort(id_vec.begin(), id_vec.end());

  // generate a string representation of the ID range(s)
  std::stringstream out;

  int i = 0;
  int start_id = id_vec[0]; // initialize with first ID
  int stop_id;
  int n_ents = id_vec.size();
  // loop over all cells in the universe
  while (i < n_ents) {

    stop_id = id_vec[i];

    // if the next ID is not in this contiguous set of IDS,
    // figure out how to write the string representing this set
    if (id_vec[i + 1] > stop_id + 1) {

      if (start_id != stop_id) {
        // there are several IDs in a row, print condensed version (i.e. 1-10,
        // 12-20)
        out << start_id << "-" << stop_id;
      } else {
        // only one ID in this contiguous block (i.e. 3, 5, 7, 9)
        out << start_id;
      }
      // insert a comma as long as we aren't in the last ID set
      if (i < n_ents - 1) {
        out << ", ";
      }

      // if we are at the end of a set, set the start ID to the first value
      // in the next set.
      start_id = id_vec[++i];
    }

    i++;
  }

  return out.str();
}

int32_t XDGUniverse::implicit_complement_idx() const
{
  xdg::MeshID ipc = xdg_ptr()->mesh_manager()->implicit_complement();
  if (ipc == xdg::ID_NONE)
    fatal_error("Implicit complement not found in XDG model");
  return cell_index(ipc);
}

bool XDGUniverse::find_cell(GeometryState& p) const
{
  // if the particle isn't in any of the other XDG
  // cells, place it in the implicit complement
  bool found = Universe::find_cell(p);
  if (!found && model::universe_map[this->id_] != model::root_universe) {
    p.lowest_coord().cell = implicit_complement_idx();
    found = true;
  }
  return found;
}

void XDGUniverse::to_hdf5(hid_t universes_group) const
{
  // Create a group for this universe.
  auto group = create_group(universes_group, fmt::format("universe {}", id_));

  // Write the geometry representation type.
  write_string(group, "geom_type", "XDG", false);

  // Write other properties of the XDG Universe
  write_string(group, "filename", filename_, false);
  write_attribute(
    group, "auto_geom_ids", static_cast<int>(adjust_geometry_ids_));
  write_attribute(
    group, "auto_mat_ids", static_cast<int>(adjust_material_ids_));

  close_group(group);
}

void XDGUniverse::assign_material(
  std::string& mat_string, std::unique_ptr<XDGCell>& c) const
{
  // material void checks
  if (mat_string == "void" || mat_string == "vacuum") {
    c->material_.push_back(MATERIAL_VOID);
    return;
  }

  bool mat_found_by_name = false;
  // attempt to find a material with a matching name
  to_lower(mat_string);
  for (const auto& m : model::materials) {
    std::string m_name = m->name();
    to_lower(m_name);
    if (mat_string == m_name) {
      // assign the material with that name
      if (!mat_found_by_name) {
        mat_found_by_name = true;
        c->material_.push_back(m->id_);
        // report error if more than one material is found
      } else {
        fatal_error(fmt::format(
          "More than one material found with name '{}'. Please ensure "
          "materials "
          "have unique names if using this property to assign materials.",
          mat_string));
      }
    }
  }

  // if no material was set using a name, assign by id
  if (!mat_found_by_name) {
    bool found_by_id = true;
    try {
      auto id = std::stoi(mat_string);
      if (model::material_map.find(id) == model::material_map.end())
        found_by_id = false;
      c->material_.emplace_back(id);
    } catch (const std::invalid_argument&) {
      found_by_id = false;
    }

    // report failure for failed int conversion or missing material
    if (!found_by_id)
      fatal_error(
        fmt::format("Material with name/ID '{}' not found for volume (cell) {}",
          mat_string, c->id_));
  }

  if (settings::verbosity >= 10) {
    const auto& m = model::materials[model::material_map.at(c->material_[0])];
    std::stringstream msg;
    msg << "XDG material " << mat_string << " was assigned";
    if (mat_found_by_name) {
      msg << " using material name: " << m->name_;
    } else {
      msg << " using material id: " << m->id_;
    }
    write_message(msg.str(), 10);
  }
}

//==============================================================================
// XDG Cell implementation
//==============================================================================

XDGCell::XDGCell(std::shared_ptr<xdg::XDG> xdg_ptr, xdg::MeshID xdg_id)
  : Cell {}, XDGGeometryObject(xdg_ptr, xdg_id)
{
  geom_type_ = GeometryType::XDG;
  // TODO: Allow XDG cells to be filled with other geometry
  fill_ = C_NONE;
};

std::pair<double, int32_t> XDGCell::distance(
  Position r, Direction u, int32_t on_surface, GeometryState* p) const
{
  // if we've changed direction or we're not on a surface,
  // reset the history and update last direction
  if (u != p->last_dir()) {
    p->last_dir() = u;
    p->xdg_prev_elements().clear();
  }
  if (on_surface == 0) {
    p->xdg_prev_elements().clear();
  }

  const auto& univ = model::universes[p->lowest_coord().universe];

  XDGUniverse* xdg_univ = static_cast<XDGUniverse*>(univ.get());
  if (!xdg_univ)
    fatal_error("XDG call made for particle in a non-XDG universe");

  // initialize to lost particle conditions
  int surf_idx = -1;
  double dist = INFINITY;


  // create the ray
  double pnt[3] = {r.x, r.y, r.z};
  double dir[3] = {u.x, u.y, u.z};
  std::pair<double, xdg::MeshID> result = xdg_ptr()->ray_fire(xdg_id(), pnt, dir, &p->xdg_prev_elements());
  if (result.second > 0) {
    surf_idx = xdg_univ->surface_index(result.second);
    dist = result.first;
  } else if (xdg_id() != xdg_ptr()->mesh_manager()->implicit_complement() || is_root_universe(xdg_univ->id_)) {
    // surface boundary conditions are ignored for projection plotting, meaning
    // that the particle may move through the graveyard (bounding) volume and
    // into the implicit complement on the other side where no intersection will
    // be found. Treating this as a lost particle is problematic when plotting.
    // Instead, the infinite distance and invalid surface index are returned.
    if (settings::run_mode == RunMode::PLOTTING)
      return {INFTY, -1};

    // the particle should be marked as lost immediately if an intersection
    // isn't found in a volume that is not the implicit complement. In the case
    // that the XDG model is the root universe of the geometry, even a missing
    // intersection in the implicit complement should trigger this condition.
    std::string material_id =
      p->material() == MATERIAL_VOID
        ? "-1 (VOID)"
        : std::to_string(model::materials[p->material()]->id());
    p->mark_as_lost(fmt::format(
      "No intersection found with XDG cell {}, filled with material {}", id_,
      material_id));
  }
  return {dist, surf_idx};
}

bool XDGCell::contains(Position r, Direction u, int32_t on_surface) const
{
  double pnt[3] = {r.x, r.y, r.z};
  xdg::Direction dir {u.x, u.y, u.z};
  return xdg_ptr()->point_in_volume(xdg_id(), pnt, &dir);
}

void XDGCell::to_hdf5_inner(hid_t group_id) const
{
  write_string(group_id, "geom_type", "XDG", false);
}

BoundingBox XDGCell::bounding_box() const
{
  xdg::BoundingBox bbox = xdg_ptr()->mesh_manager()->volume_bounding_box(xdg_id());
  return {bbox[0], bbox[1], bbox[2], bbox[3], bbox[4], bbox[5]};
}

//==============================================================================
// XDGSurface implementation
//==============================================================================

XDGSurface::XDGSurface(std::shared_ptr<xdg::XDG> xdg_ptr, int32_t xdg_id)
  : Surface {}, XDGGeometryObject(xdg_ptr, xdg_id)
{
  geom_type_ = GeometryType::XDG;
}

double XDGSurface::evaluate(Position r) const
{
  return 0.0;
}

double XDGSurface::distance(Position r, Direction u, bool coincident) const
{
  double pnt[3] = {r.x, r.y, r.z};
  double dir[3] = {u.x, u.y, u.z};
  std::pair<double, xdg::MeshID> result = xdg_ptr()->ray_fire_surface(xdg_id(), pnt, dir);
  return result.first < 0.0 ? INFTY : result.first;
}

Direction XDGSurface::normal(Position r) const
{
  xdg::Direction normal = xdg_ptr()->surface_normal(xdg_id(), {r.x, r.y, r.z});
  return {normal[0], normal[1], normal[2]};
}

Direction XDGSurface::reflect(Position r, Direction u, GeometryState* p) const
{
  Expects(p);
  p->xdg_prev_elements() = {p->xdg_prev_elements().back()};
  double pnt[3] = {r.x, r.y, r.z};
  xdg::Direction normal = xdg_ptr()->surface_normal(xdg_id(), pnt, &p->xdg_prev_elements());
  p->last_dir() = u.reflect({normal[0], normal[1], normal[2]});
  return p->last_dir();
}

//==============================================================================
// Non-member functions
//==============================================================================

void read_xdg_universes(pugi::xml_node node)
{
  for (pugi::xml_node dag_node : node.children("xdg_universe")) {
    model::universes.push_back(std::make_unique<XDGUniverse>(dag_node));
    model::universe_map[model::universes.back()->id_] =
      model::universes.size() - 1;
  }
}

void check_xdg_root_univ()
{
  return;
  //// TODO: Find an alternative method for this check
  // const auto& ru = model::universes[model::root_universe];
  // if (ru->geom_type() == GeometryType::XDG) {
  //   // if the root universe contains XDG geometry, warn the user
  //   // if it does not contain a graveyard volume
  //   auto dag_univ = dynamic_cast<XDGUniverse*>(ru.get());
  //   if (dag_univ && !dag_univ->has_graveyard()) {
  //     warning(
  //       "No graveyard volume found in the XDG model. "
  //       "This may result in lost particles and rapid simulation failure.");
  //   }
  // }
}

int32_t xdg_next_cell(int32_t surf, int32_t curr_cell, int32_t univ)
{
  auto surfp = dynamic_cast<XDGSurface*>(model::surfaces[surf - 1].get());
  auto cellp = dynamic_cast<XDGCell*>(model::cells[curr_cell].get());
  auto univp = static_cast<XDGUniverse*>(model::universes[univ].get());

  xdg::MeshID next_volume = univp->xdg_ptr()->next_volume(cellp->xdg_id(), surfp->xdg_id());
  return univp->cell_index(next_volume);
}

} // namespace openmc

#else

namespace openmc {

void read_xdg_universes(pugi::xml_node node)
{
  // TODO : updated to xdg universe
  if (check_for_node(node, "xdg_universe")) {
    fatal_error("XDG Universes are present but OpenMC was not configured "
                "with XDG");
  }
};

void check_xdg_root_univ() {};

int32_t next_cell(int32_t surf, int32_t curr_cell, int32_t univ);

} // namespace openmc

#endif // XDG
