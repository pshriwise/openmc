#include "openmc/universe.h"

#include <set>

#include "openmc/hdf5_interface.h"
#include "openmc/material.h"
#include "openmc/particle.h"
#include "openmc/string_utils.h"

namespace openmc {

namespace model {

std::unordered_map<int32_t, int32_t> universe_map;
vector<unique_ptr<Universe>> universes;

} // namespace model

//==============================================================================
// Universe implementation
//==============================================================================

void Universe::to_hdf5(hid_t universes_group) const
{
  // Create a group for this universe.
  auto group = create_group(universes_group, fmt::format("universe {}", id_));

  // Write the geometry representation type.
  write_string(group, "geom_type", "csg", false);

  // Write the contained cells.
  if (cells_.size() > 0) {
    vector<int32_t> cell_ids;
    for (auto i_cell : cells_)
      cell_ids.push_back(model::cells[i_cell]->id_);
    write_dataset(group, "cells", cell_ids);
  }

  close_group(group);
}

bool Universe::find_cell(GeometryState& p) const
{
  const auto& cells {
    !partitioner_ ? cells_ : partitioner_->get_cells(p.r_local(), p.u_local())};

  Position r {p.r_local()};
  Position u {p.u_local()};
  auto surf = p.surface();
  int32_t i_univ = p.lowest_coord().universe;

  for (auto i_cell : cells) {
    if (model::cells[i_cell]->universe_ != i_univ)
      continue;
    // Check if this cell contains the particle
    if (model::cells[i_cell]->contains(r, u, surf)) {
      p.lowest_coord().cell = i_cell;
      return true;
    }
  }
  return false;
}

BoundingBox Universe::bounding_box() const
{
  BoundingBox bbox = {INFTY, -INFTY, INFTY, -INFTY, INFTY, -INFTY};
  if (cells_.size() == 0) {
    return {};
  } else {
    for (const auto& cell : cells_) {
      auto& c = model::cells[cell];
      bbox |= c->bounding_box();
    }
  }
  return bbox;
}

//==============================================================================
// MeshUniverse implementation
//==============================================================================

#ifdef LIBMESH

MeshUniverse::MeshUniverse(pugi::xml_node node)
{
  geom_type_ = GeometryType::MESH;

  if (check_for_node(node, "id")) {
    id_ = std::stoi(get_node_value(node, "id"));
  } else {
    fatal_error("Must specify the ID of the Mesh Universe");
  }

  if (check_for_node(node, "name")) {
    name_ = std::stoi(get_node_value(node, "name"));
  }

  if (check_for_node(node, "mesh")) {
    int32_t mesh_id = std::stoi(get_node_value(node, "mesh"));
    // Ensure the speicifed mesh is present
    if (model::mesh_map.find(mesh_id) == model::mesh_map.end()) {
      fatal_error(fmt::format("Mesh {} could not be found", mesh_id));
    }
    mesh_ = model::mesh_map[mesh_id];
  } else {
    fatal_error(fmt::format("No mesh specified on mesh universe {}", id_));
  }

  create_cells(node);
  set_boundary_conditions();
}

void MeshUniverse::create_unstructured_mesh_cells()
{
  // get the openmc::LibMesh mesh and libmesh mesh
  const auto& mesh = model::meshes[mesh_];
  LibMesh* mesh_ptr = dynamic_cast<LibMesh*>(mesh.get());
  const auto& lmesh = mesh_ptr->libmesh_mesh();

  std::set<libMesh::subdomain_id_type> subdomain_ids;
  lmesh->subdomain_ids(subdomain_ids);

  int mat_elements {0};

  for (const auto& subdomain_id : subdomain_ids) {
    std::string subdomain_name = lmesh->subdomain_name(subdomain_id);
    std::cout << "Subdomain name: " << subdomain_name << std::endl;
    std::string lower_name = subdomain_name;
    to_lower(lower_name);

    if (subdomain_name == "")
      continue;
    // get the elements for this subdomain id
    auto subdomain_elements =
      lmesh->active_subdomain_elements_ptr_range(subdomain_id);

    // determine the material ID to assign for this block
    int32_t subdomain_mat_id = MATERIAL_VOID;
    int32_t mat_by_name_idx = get_material_by_name(subdomain_name);
    if (lower_name == "vacuum") {
      subdomain_mat_id = MATERIAL_VOID;
    } else if (mat_by_name_idx != -1) {
      subdomain_mat_id = model::materials[mat_by_name_idx]->id();
      write_message(
        fmt::format("Assigning subdomain {} by name: {} with material ID {}.",
          subdomain_id, subdomain_name, subdomain_mat_id),
        10);
    } else if (model::material_map.find(subdomain_id) !=
               model::material_map.end()) {
      // if no matching name is found, set cell materials by domain id
      write_message(fmt::format("Assigning subdomain {} by ID.", subdomain_id,
                      subdomain_name),
        10);
      subdomain_mat_id = subdomain_id;
    } // else {
    //   fatal_error(fmt::format("Could not find material by name ({}) or ID
    //   ({}).", subdomain_name, subdomain_id));
    // }

    int n_subdomain_elems {0};
    // set each cell in the block with the chosen material ID
    for (const auto& elem_ptr : subdomain_elements) {
      int bin = mesh_ptr->get_bin_from_element(elem_ptr);
      // replace temporary void material assignment
      model::cells[cells_[bin]]->material_[0] = subdomain_mat_id;
      mat_elements++;
      n_subdomain_elems++;
    }
    std::cout << fmt::format("{} elements\n", n_subdomain_elems);
  }

  std::cout << fmt::format(
    "Assigned materials to {} elements.\n", mat_elements);
}

bool MeshUniverse::is_reflecting_face(int elem_idx, int face_idx) const
{
  const auto* lmesh = dynamic_cast<LibMesh*>(model::meshes[mesh()].get());
  const auto* elem_ptr = lmesh->get_element_ptr_from_bin(elem_idx);

  std::vector<libMesh::boundary_id_type> boundary_ids;
  lmesh->mesh_ptr()->boundary_info->boundary_ids(
    elem_ptr, face_idx, boundary_ids);

  return false;
}

void MeshUniverse::create_cells(pugi::xml_node node)
{
  vector<std::string> cell_fills;
  if (check_for_node(node, "fills")) {
    cell_fills = get_node_array<std::string>(node, "fills");
  }

  bool structured_mesh =
    model::meshes[mesh_]->structure() == MeshStructure::STRUCTURED;

  int n_bins = model::meshes[mesh_]->n_bins();
  if (cell_fills.size() != 1 && cell_fills.size() != n_bins &&
      structured_mesh) {
    fatal_error(
      fmt::format("Invalid number of cell fills provided for structured mesh "
                  "universe {}. Must be 1 or {}",
        id_, n_bins));
  }

  cells_.reserve(n_bins);
  // find the available cell id
  int32_t next_cell_id {-1};
  for (const auto& c : model::cells) {
    next_cell_id = std::max(next_cell_id, c->id_);
  }
  next_cell_id++;

  // create cells to fill the mesh elements
  int32_t fill = MATERIAL_VOID;
  for (int i = 0; i < n_bins; i++) {
    // if more than one cell fill is provided, assume that each mesh
    // element has its own fill

    if (structured_mesh) {
      if (cell_fills.size() > 1)
        fill = std::stoi(cell_fills[0]);
      else
        fill = std::stoi(cell_fills[i]);
      // check that this fill is in the material array
      if (model::material_map.find(fill) == model::material_map.end()) {
        fatal_error(
          fmt::format("Material {} not found for MeshUniverse {}", fill, id_));
      }
    }

    // create a new mesh cell
    model::cells.push_back(std::make_unique<MeshCell>(mesh_, i));
    const auto& cell = model::cells.back();

    cell->type_ = Fill::MATERIAL;
    cell->material_.push_back(fill);
    cell->id_ = next_cell_id++;
    // set universe ID, this will be updated later by another loop in
    // geometry_aux
    cell->universe_ = id_;
    model::cell_map[cell->id_] = model::cells.size() - 1;
    cells_[i] = model::cells.size() - 1;
  }

  // special function if mesh is unstructured
  if (model::meshes[mesh_]->structure() == MeshStructure::UNSTRUCTURED) {
    create_unstructured_mesh_cells();
  }

  // create a cell for the exterior of the mesh
  int32_t outer_material = MATERIAL_VOID;
  if (check_for_node(node, "outer")) {
    // get id of outer material
    outer_material = std::stoi(get_node_value(node, "outer"));
  }
  // negative one indicates that this cell is the exterior of the mesh
  model::cells.push_back(std::make_unique<MeshCell>(mesh_, -1));
  const auto& cell = model::cells.back();

  cell->type_ = Fill::MATERIAL;
  cell->material_.push_back(outer_material);
  cell->id_ = next_cell_id;
  cell->universe_ = id_;

  model::cell_map[cell->id_] = model::cells.size() - 1;
  outer() = model::cell_map[cell->id_];
}

bool MeshUniverse::find_cell(Particle& p) const
{
  Position r {p.r_local()};
  const auto& mesh = model::meshes[mesh_];
  bool in_mesh;
  int mesh_bin = mesh->get_bin(r);

  if (mesh_bin == C_NONE) {
    if (outer() == C_NONE)
      return false;
    p.coord(p.n_coord() - 1).mesh_cell_index() = mesh_bin;
    p.coord(p.n_coord() - 1).cell = outer();
    return true;
  }

  p.coord(p.n_coord() - 1).mesh_cell_index() = mesh_bin;
  p.coord(p.n_coord() - 1).cell = cells_[mesh_bin];
  return true;
}

void MeshUniverse::set_boundary_conditions()
{
  const auto* lmesh = libmesh_ptr()->mesh_ptr();
  const auto& boundary_info = lmesh->boundary_info;

  // loop over all boundary sets and print name
  std::cout << "Number of boundary sets: " << boundary_info->n_boundary_ids()
            << std::endl;
  for (auto sideset_info : boundary_info->set_sideset_name_map()) {
    std::cout << "Sideset " << sideset_info.first
              << " name: " << sideset_info.second << std::endl;
  }
}

void MeshUniverse::next_cell(Particle& p) const
{
  auto& coord = p.coord(p.n_coord() - 1);
  const auto mesh = dynamic_cast<LibMesh*>(model::meshes[mesh_].get());

  int32_t next_mesh_idx = p.boundary().mesh_translation(0);
  int32_t next_cell_idx {C_NONE};
  if (mesh->bin_is_valid(next_mesh_idx)) {

    if (p.coord(p.n_coord() - 1).mesh_cell_index() == C_NONE) {
      write_message(
        fmt::format(
          "\tParticle {} moving into the mesh. \n\tPosition: {} {} {}", p.id(),
          p.r()[0], p.r()[1], p.r()[2]),
        10);
    }
    write_message(
      fmt::format("\nMoving into mesh cell: {} {} {}",
        p.boundary().mesh_translation(0), p.boundary().mesh_translation(1),
        p.boundary().mesh_translation(2)),
      10);

    next_cell_idx = cells_[next_mesh_idx];
  } else {
    write_message(
      fmt::format("\tParticle {} moving out of the mesh. \n\tPosition: {} {} "
                  "{} \n\tDirection: {} {} {}",
        p.id(), p.r()[0], p.r()[1], p.r()[2], p.u()[0], p.u()[1], p.u()[2]),
      10);
    p.wgt() = 0.0;
    next_mesh_idx = C_NONE;
    next_cell_idx = outer_;
  }

  // reset the lattice_translation for the boundary crossing
  p.boundary().mesh_translation() = {0, 0, 0};
  // update material and temperature
  p.material_last() = p.material();
  p.sqrtkT_last() = p.sqrtkT();
  // set previous bin
  coord.mesh_cell_index() = next_mesh_idx;
  coord.cell = next_cell_idx;
  coord.mesh_index() = p.boundary().mesh_translation();
  const auto& cell = model::cells.at(next_cell_idx);
  Cell* cell_ptr = model::cells.at(next_cell_idx).get();
  // TODO: Support multiple cell instances
  p.cell_instance() = 0;
  p.material() = cell->material_[0];
  p.sqrtkT() = cell->sqrtkT_[0];
}

#endif

//==============================================================================
// UniversePartitioner implementation
//==============================================================================

UniversePartitioner::UniversePartitioner(const Universe& univ)
{
  // Define an ordered set of surface indices that point to z-planes.  Use a
  // functor to to order the set by the z0_ values of the corresponding planes.
  struct compare_surfs {
    bool operator()(const int32_t& i_surf, const int32_t& j_surf) const
    {
      const auto* surf = model::surfaces[i_surf].get();
      const auto* zplane = dynamic_cast<const SurfaceZPlane*>(surf);
      double zi = zplane->z0_;
      surf = model::surfaces[j_surf].get();
      zplane = dynamic_cast<const SurfaceZPlane*>(surf);
      double zj = zplane->z0_;
      return zi < zj;
    }
  };
  std::set<int32_t, compare_surfs> surf_set;

  // Find all of the z-planes in this universe.  A set is used here for the
  // O(log(n)) insertions that will ensure entries are not repeated.
  for (auto i_cell : univ.cells_) {
    for (auto token : model::cells[i_cell]->surfaces()) {
      auto i_surf = std::abs(token) - 1;
      const auto* surf = model::surfaces[i_surf].get();
      if (const auto* zplane = dynamic_cast<const SurfaceZPlane*>(surf))
        surf_set.insert(i_surf);
    }
  }

  // Populate the surfs_ vector from the ordered set.
  surfs_.insert(surfs_.begin(), surf_set.begin(), surf_set.end());

  // Populate the partition lists.
  partitions_.resize(surfs_.size() + 1);
  for (auto i_cell : univ.cells_) {
    // It is difficult to determine the bounds of a complex cell, so add complex
    // cells to all partitions.
    if (!model::cells[i_cell]->is_simple()) {
      for (auto& p : partitions_)
        p.push_back(i_cell);
      continue;
    }

    // Find the tokens for bounding z-planes.
    int32_t lower_token = 0, upper_token = 0;
    double min_z, max_z;
    for (auto token : model::cells[i_cell]->surfaces()) {
      const auto* surf = model::surfaces[std::abs(token) - 1].get();
      if (const auto* zplane = dynamic_cast<const SurfaceZPlane*>(surf)) {
        if (lower_token == 0 || zplane->z0_ < min_z) {
          lower_token = token;
          min_z = zplane->z0_;
        }
        if (upper_token == 0 || zplane->z0_ > max_z) {
          upper_token = token;
          max_z = zplane->z0_;
        }
      }
    }

    // If there are no bounding z-planes, add this cell to all partitions.
    if (lower_token == 0) {
      for (auto& p : partitions_)
        p.push_back(i_cell);
      continue;
    }

    // Find the first partition this cell lies in.  If the lower_token indicates
    // a negative halfspace, then the cell is unbounded in the lower direction
    // and it lies in the first partition onward.  Otherwise, it is bounded by
    // the positive halfspace given by the lower_token.
    int first_partition = 0;
    if (lower_token > 0) {
      for (int i = 0; i < surfs_.size(); ++i) {
        if (lower_token == surfs_[i] + 1) {
          first_partition = i + 1;
          break;
        }
      }
    }

    // Find the last partition this cell lies in.  The logic is analogous to the
    // logic for first_partition.
    int last_partition = surfs_.size();
    if (upper_token < 0) {
      for (int i = first_partition; i < surfs_.size(); ++i) {
        if (upper_token == -(surfs_[i] + 1)) {
          last_partition = i;
          break;
        }
      }
    }

    // Add the cell to all relevant partitions.
    for (int i = first_partition; i <= last_partition; ++i) {
      partitions_[i].push_back(i_cell);
    }
  }
}

const vector<int32_t>& UniversePartitioner::get_cells(
  Position r, Direction u) const
{
  // Perform a binary search for the partition containing the given coordinates.
  int left = 0;
  int middle = (surfs_.size() - 1) / 2;
  int right = surfs_.size() - 1;
  while (true) {
    // Check the sense of the coordinates for the current surface.
    const auto& surf = *model::surfaces[surfs_[middle]];
    if (surf.sense(r, u)) {
      // The coordinates lie in the positive halfspace.  Recurse if there are
      // more surfaces to check.  Otherwise, return the cells on the positive
      // side of this surface.
      int right_leaf = right - (right - middle) / 2;
      if (right_leaf != middle) {
        left = middle + 1;
        middle = right_leaf;
      } else {
        return partitions_[middle + 1];
      }

    } else {
      // The coordinates lie in the negative halfspace.  Recurse if there are
      // more surfaces to check.  Otherwise, return the cells on the negative
      // side of this surface.
      int left_leaf = left + (middle - left) / 2;
      if (left_leaf != middle) {
        right = middle - 1;
        middle = left_leaf;
      } else {
        return partitions_[middle];
      }
    }
  }
}

void read_mesh_universes(pugi::xml_node node)
{
#ifdef LIBMESH
  int n_univs = 0;
  for (pugi::xml_node mesh_univ_node : node.children("mesh_universe")) {
    model::universes.push_back(std::make_unique<MeshUniverse>(mesh_univ_node));
    model::universe_map[model::universes.back()->id_] = n_univs++;
  }
#endif
}

} // namespace openmc
