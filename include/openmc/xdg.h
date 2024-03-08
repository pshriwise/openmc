#ifndef OPENMC_XDG_H
#define OPENMC_XDG_H

namespace openmc {
extern "C" const bool XDG_ENABLED;
}

// always include the XML interface header
#include "openmc/xml_interface.h"

//==============================================================================
// Functions that are always defined
//==============================================================================

namespace openmc {

void read_dagmc_universes(pugi::xml_node node);
void check_dagmc_root_univ();

} // namespace openmc

#ifdef OPENMC_XDG

#include "xdg/xdg.h"

#include "openmc/cell.h"
#include "openmc/particle.h"
#include "openmc/position.h"
#include "openmc/surface.h"

namespace openmc {

class XDGGeometryObject {
public:
  XDGGeometryObject(std::shared_ptr<xdg::XDG> xdg_ptr, xdg::MeshID xdg_id) :
    xdg_ptr_{xdg_ptr}, xdg_id_{xdg_id} {}

  // Accessor methods
  xdg::MeshID xdg_id() const { return xdg_id_; }
  const xdg::XDG* xdg_ptr() const { return xdg_ptr_.get(); }

private:
  std::shared_ptr<xdg::XDG> xdg_ptr_;      //!< Pointer to XDG instance
  xdg::MeshID xdg_id_;                     //!< XDG ID
};

class XDGSurface : public XDGGeometryObject, public Surface {
public:
  XDGSurface(std::shared_ptr<xdg::XDG> dag_ptr, xdg::MeshID xdg_id);

  double evaluate(Position r) const override;
  double distance(Position r, Direction u, bool coincident) const override;

  Direction normal(Position r) const override;
  Direction reflect(Position r, Direction u, GeometryState* p) const override;

  inline void to_hdf5_inner(hid_t group_id) const override {};
};

class XDGCell : public XDGGeometryObject, public Cell {
public:
  XDGCell(std::shared_ptr<xdg::XDG> xdg_ptr, xdg::MeshID xdg_id);

  bool contains(Position r, Direction u, int32_t on_surface) const override;

  std::pair<double, int32_t> distance(Position r, Direction u,
    int32_t on_surface, GeometryState* p) const override;

  BoundingBox bounding_box() const override;

  void to_hdf5_inner(hid_t group_id) const override;
};

class XDGUniverse : public Universe {

public:
  explicit XDGUniverse(pugi::xml_node node);

  //! Create a new XDG universe
  //! \param[in] filename Name of the XDG file
  //! \param[in] auto_geom_ids Whether or not to automatically assign cell and
  //! surface IDs
  //! \param[in] auto_mat_ids Whether or not to automatically assign
  //! material IDs
  explicit XDGUniverse(const std::string& filename, bool auto_geom_ids = false,
    bool auto_mat_ids = false);

  //! Alternative XDG universe constructor for external XDG instance
  explicit XDGUniverse(std::shared_ptr<xdg::XDG> external_xdg_ptr,
    const std::string& filename = "", bool auto_geom_ids = false,
    bool auto_mat_ids = false);

  //! Initialize the XDG accel. data structures, indices, material
  //! assignments, etc.
  void initialize();

  //! Returns the index to the implicit complement's index in OpenMC for this
  //! XDG universe
  int32_t implicit_complement_idx() const;

  //! Assign a material to a cell based
  //! \param[in] mat_string The XDG material assignment string
  //! \param[in] c The OpenMC cell to which the material is assigned
  void assign_material(
    std::string& mat_string, std::unique_ptr<XDGCell>& c) const;

  //! Return the index into the model cells vector for a given XDG volume
  //! handle in the universe
  //! \param[in] vol MOAB handle to the XDG volume set
  int32_t cell_index(xdg::MeshID voume) const;

  //! Return the index into the model surfaces vector for a given XDG surface
  //! handle in the universe
  //! \param[in] surf MOAB handle to the XDG surface set
  int32_t surface_index(xdg::MeshID surface) const;

  //! Generate a string representing the ranges of IDs present in the XDG
  //! model. Contiguous chunks of IDs are represented as a range (i.e. 1-10). If
  //! there is a single ID a chunk, it will be represented as a single number
  //! (i.e. 2, 4, 6, 8).
  //! \param[in] dim Dimension of the entities
  //! \return A string of the ID ranges for entities of dimension \p dim
  std::string xdg_ids_for_dim(int dim) const;

  bool find_cell(GeometryState& p) const override;

  void to_hdf5(hid_t universes_group) const override;

  // Data Members
  std::shared_ptr<xdg::XDG>
    xdg_instance_;        //!< XDG Instance for this universe
  int32_t cell_idx_offset_; //!< An offset to the start of the cells in this
                            //!< universe in OpenMC's cell vector
  int32_t surf_idx_offset_; //!< An offset to the start of the surfaces in this
                            //!< universe in OpenMC's surface vector

  // Accessors
  xdg::XDG* xdg_ptr() const { return xdg_instance_.get(); }

private:
  void set_id();        //!< Deduce the universe id from model::universes
  void init_xdg();    //!< Create and initialise XDG pointer
  void init_metadata(); //!< Create and initialise dagmcMetaData pointer
  void init_geometry(); //!< Create cells and surfaces from XDG entities

  std::string
    filename_; //!< Name of the XDG file used to create this universe

  bool adjust_geometry_ids_; //!< Indicates whether or not to automatically
                             //!< generate new cell and surface IDs for the
                             //!< universe
  bool adjust_material_ids_; //!< Indicates whether or not to automatically
                             //!< generate new material IDs for the universe

  // mappings from XDG IDs to OpenMC surface and cell indices
  std::unordered_map<xdg::MeshID, int32_t> surface_index_map_;
  std::unordered_map<xdg::MeshID, int32_t> cell_index_map_;
};

//==============================================================================
// Non-member functions
//==============================================================================

int32_t xdg_next_cell(int32_t surf, int32_t curr_cell, int32_t univ);

} // namespace openmc

#endif // XDG

#endif // OPENMC_XDG_H
