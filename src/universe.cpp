#include "openmc/universe.h"

#include <fmt/core.h>
#include <set>

#include "openmc/hdf5_interface.h"
#include "openmc/material.h"
#include "openmc/mesh.h"
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

bool Universe::find_cell(Particle& p) const
{
  const auto& cells {
    !partitioner_ ? cells_ : partitioner_->get_cells(p.r_local(), p.u_local())};

  for (auto it = cells.begin(); it != cells.end(); it++) {
    int32_t i_cell = *it;
    int32_t i_univ = p.coord(p.n_coord() - 1).universe;
    if (model::cells[i_cell]->universe_ != i_univ)
      continue;

    // Check if this cell contains the particle;
    Position r {p.r_local()};
    Direction u {p.u_local()};
    auto surf = p.surface();
    if (model::cells[i_cell]->contains(r, u, surf)) {
      p.coord(p.n_coord() - 1).cell = i_cell;
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

  if (check_for_node(node, "fills")) {
    auto fill_strs = get_node_array<std::string>(node, "fills");
    create_cells(fill_strs);
  }

  if (check_for_node(node, "outer")) {
    outer() = std::stoi(get_node_value(node, "outer"));
  }
}

void MeshUniverse::create_cells(const vector<std::string>& cell_fills)
{
  int n_bins = model::meshes[mesh_]->n_bins();
  if (cell_fills.size() != 1 && cell_fills.size() != n_bins) {
    fatal_error(fmt::format("Invalid number of cell fills provided for mesh "
                            "universe {}. Must be 1 or {}",
      id_, n_bins));
  }

  cells_.reserve(n_bins);
  // find the available cell id
  int32_t next_cell_id;
  for (const auto& c : model::cells) {
    next_cell_id = std::max(next_cell_id, c->id_);
  }
  next_cell_id++;

  // create cells to fill the mesh elements
  // TODO: extend beyond material fills
  int32_t fill = std::stoi(cell_fills[0]);
  for (int i = 0; i < n_bins; i++) {
    // if more than one cell fill is provided, assume that each mesh
    // element has its own fill
    if (cell_fills.size() != 1)
      fill = std::stoi(cell_fills[i]);
    // check that this fill is in the material array
    if (model::material_map.find(fill) == model::material_map.end()) {
      fatal_error(
        fmt::format("Material {} not found for MeshUniverse {}", fill, id_));
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
}

bool MeshUniverse::find_cell(Particle& p) const
{
  Position r {p.r_local()};

  int mesh_bin = model::meshes[mesh_]->get_bin(r);
  if (mesh_bin == -1) {
    if (outer() == C_NONE)
      return false;
    p.coord(p.n_coord() - 1).universe = outer();
    return model::universes[outer()]->find_cell(p);
  }

  p.coord(p.n_coord() - 1).mesh_cell_index() = mesh_bin;
  p.coord(p.n_coord() - 1).cell = cells_[mesh_bin];
  return true;
}

void MeshUniverse::next_cell(Particle& p) const
{
  auto& coord = p.coord(p.n_coord() - 1);
  const auto mesh = dynamic_cast<StructuredMesh*>(model::meshes[mesh_].get());

  // the particle's surface attributte should contain the flat index
  // of the mesh cell it will enter next
  int32_t next_mesh_idx =
    mesh->get_bin_from_indices(p.boundary().lattice_translation);
  int32_t next_cell_idx = cells_[next_mesh_idx];

  p.boundary().lattice_translation[0] = 0;
  p.boundary().lattice_translation[1] = 0;
  p.boundary().lattice_translation[2] = 0;
  // update material and temperature
  p.material_last() = p.material();
  p.sqrtkT_last() = p.sqrtkT();
  // set new cell value
  p.coord(p.n_coord() - 1).mesh_cell_index() = next_mesh_idx;
  p.coord(p.n_coord() - 1).cell = next_cell_idx;
  p.coord(p.n_coord() - 1).lattice_i = p.boundary().lattice_translation;
  const auto& cell = model::cells.at(next_cell_idx);
  // TODO: Support multiple cell instances
  p.cell_instance() = 0;
  p.material() = cell->material_[0];
  p.sqrtkT() = cell->sqrtkT_[0];
  return;
}

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
    for (auto token : model::cells[i_cell]->rpn_) {
      if (token < OP_UNION) {
        auto i_surf = std::abs(token) - 1;
        const auto* surf = model::surfaces[i_surf].get();
        if (const auto* zplane = dynamic_cast<const SurfaceZPlane*>(surf))
          surf_set.insert(i_surf);
      }
    }
  }

  // Populate the surfs_ vector from the ordered set.
  surfs_.insert(surfs_.begin(), surf_set.begin(), surf_set.end());

  // Populate the partition lists.
  partitions_.resize(surfs_.size() + 1);
  for (auto i_cell : univ.cells_) {
    // It is difficult to determine the bounds of a complex cell, so add complex
    // cells to all partitions.
    if (!model::cells[i_cell]->simple_) {
      for (auto& p : partitions_)
        p.push_back(i_cell);
      continue;
    }

    // Find the tokens for bounding z-planes.
    int32_t lower_token = 0, upper_token = 0;
    double min_z, max_z;
    for (auto token : model::cells[i_cell]->rpn_) {
      if (token < OP_UNION) {
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
  vector<int> mesh_universe_ids;
  for (pugi::xml_node mesh_univ_node : node.children("mesh_universe")) {
    model::universes.push_back(std::make_unique<MeshUniverse>(mesh_univ_node));
    mesh_universe_ids.push_back(model::universes.back()->id_);
    model::universe_map[model::universes.back()->id_] =
      model::universes.size() - 1;
  }
}

} // namespace openmc