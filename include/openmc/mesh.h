//! \file mesh.h
//! \brief Mesh types used for tallies, Shannon entropy, CMFD, etc.

#ifndef OPENMC_MESH_H
#define OPENMC_MESH_H

#include <memory> // for unique_ptr
#include <vector>
#include <unordered_map>

#include "hdf5.h"
#include "pugixml.hpp"
#include "xtensor/xtensor.hpp"

#include "openmc/particle.h"
#include "openmc/position.h"

#ifdef DAGMC
#include "moab/Core.hpp"
#include "moab/AdaptiveKDTree.hpp"
#include "moab/Matrix3.hpp"
#include "moab/GeomUtil.hpp"
#endif

#ifdef LIBMESH
#include "libmesh/libmesh.h"
#include "libmesh/elem.h"
#include "libmesh/equation_systems.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/explicit_system.h"
#include "libmesh/dof_map.h"
#include "libmesh/mesh.h"
#include "libmesh/point.h"
#endif

namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

class Mesh;

namespace model {

extern std::vector<std::unique_ptr<Mesh>> meshes;
extern std::unordered_map<int32_t, int32_t> mesh_map;

} // namespace model


class Mesh
{
public:
  // Constructors and destructor
  Mesh() = default;
  Mesh(pugi::xml_node node);
  virtual ~Mesh() = default;

  // Methods

  //! Determine which bins were crossed by a particle
  //
  //! \param[in] p Particle to check
  //! \param[out] bins Bins that were crossed
  //! \param[out] lengths Fraction of tracklength in each bin
  virtual void bins_crossed(const Particle* p, std::vector<int>& bins,
                            std::vector<double>& lengths) const = 0;

  //! Determine which surface bins were crossed by a particle
  //
  //! \param[in] p Particle to check
  //! \param[out] bins Surface bins that were crossed
  virtual void
  surface_bins_crossed(const Particle* p, std::vector<int>& bins) const = 0;

  //! Get bin at a given position in space
  //
  //! \param[in] r Position to get bin for
  //! \return Mesh bin
  virtual int get_bin(Position r) const = 0;

  //! Get bin given mesh indices
  //
  //! \param[in] Array of mesh indices
  //! \return Mesh bin
  virtual int get_bin_from_indices(const int* ijk) const = 0;

  //! Get mesh indices given a position
  //
  //! \param[in] r Position to get indices for
  //! \param[out] ijk Array of mesh indices
  //! \param[out] in_mesh Whether position is in mesh
  virtual void get_indices(Position r, int* ijk, bool* in_mesh) const = 0;

  //! Get mesh indices corresponding to a mesh bin
  //
  //! \param[in] bin Mesh bin
  //! \param[out] ijk Mesh indices
  virtual void get_indices_from_bin(int bin, int* ijk) const = 0;

  //! Get the number of mesh cells.
  virtual int n_bins() const = 0;

  //! Get the number of mesh cell surfaces.
  virtual int n_surface_bins() const = 0;

  //! Find the mesh lines that intersect an axis-aligned slice plot
  //
  //! \param[in] plot_ll The lower-left coordinates of the slice plot.
  //! \param[in] plot_ur The upper-right coordinates of the slice plot.
  //! \return A pair of vectors indicating where the mesh lines lie along each
  //!   of the plot's axes.  For example an xy-slice plot will get back a vector
  //!   of x-coordinates and another of y-coordinates.  These vectors may be
  //!   empty for low-dimensional meshes.
  virtual std::pair<std::vector<double>, std::vector<double>>
  plot(Position plot_ll, Position plot_ur) const = 0;

  //! Write mesh data to an HDF5 group
  //
  //! \param[in] group HDF5 group
  virtual void to_hdf5(hid_t group) const = 0;

  // Data members

  int id_ {-1};  //!< User-specified ID
  int n_dimension_; //!< Number of dimensions
  xt::xtensor<double, 1> lower_left_; //!< Lower-left coordinates of mesh
  xt::xtensor<double, 1> upper_right_; //!< Upper-right coordinates of mesh
};

//==============================================================================
//! Tessellation of n-dimensional Euclidean space by congruent squares or cubes
//==============================================================================

class RegularMesh : public Mesh
{
public:
  // Constructors
  RegularMesh() = default;
  RegularMesh(pugi::xml_node node);

  // Overriden methods

  void bins_crossed(const Particle* p, std::vector<int>& bins,
                    std::vector<double>& lengths) const override;

  void surface_bins_crossed(const Particle* p, std::vector<int>& bins)
  const override;

  int get_bin(Position r) const override;

  int get_bin_from_indices(const int* ijk) const override;

  void get_indices(Position r, int* ijk, bool* in_mesh) const override;

  void get_indices_from_bin(int bin, int* ijk) const override;

  int n_bins() const override;

  int n_surface_bins() const override;

  std::pair<std::vector<double>, std::vector<double>>
  plot(Position plot_ll, Position plot_ur) const override;

  void to_hdf5(hid_t group) const override;

  // New methods

  //! Check where a line segment intersects the mesh and if it intersects at all
  //
  //! \param[in,out] r0 In: starting position, out: intersection point
  //! \param[in] r1 Ending position
  //! \param[out] ijk Indices of the mesh bin containing the intersection point
  //! \return Whether the line segment connecting r0 and r1 intersects mesh
  bool intersects(Position& r0, Position r1, int* ijk) const;


  //! Count number of bank sites in each mesh bin / energy bin
  //
  //! \param[in] bank Array of bank sites
  //! \param[out] Whether any bank sites are outside the mesh
  //! \return Array indicating number of sites in each mesh/energy bin
  xt::xtensor<double, 1> count_sites(const std::vector<Particle::Bank>& bank,
    bool* outside) const;

  int num_bins() const;

  // Data members

  double volume_frac_; //!< Volume fraction of each mesh element
  xt::xtensor<int, 1> shape_; //!< Number of mesh elements in each dimension
  xt::xtensor<double, 1> width_; //!< Width of each mesh element

private:
  bool intersects_1d(Position& r0, Position r1, int* ijk) const;
  bool intersects_2d(Position& r0, Position r1, int* ijk) const;
  bool intersects_3d(Position& r0, Position r1, int* ijk) const;
};


class RectilinearMesh : public Mesh
{
public:
  // Constructors
  RectilinearMesh(pugi::xml_node node);

  // Overriden methods

  void bins_crossed(const Particle* p, std::vector<int>& bins,
                    std::vector<double>& lengths) const override;

  void surface_bins_crossed(const Particle* p, std::vector<int>& bins)
  const override;

  int get_bin(Position r) const override;

  int get_bin_from_indices(const int* ijk) const override;

  void get_indices(Position r, int* ijk, bool* in_mesh) const override;

  void get_indices_from_bin(int bin, int* ijk) const override;

  int n_bins() const override;

  int n_surface_bins() const override;

  std::pair<std::vector<double>, std::vector<double>>
  plot(Position plot_ll, Position plot_ur) const override;

  void to_hdf5(hid_t group) const override;

  // New methods

  //! Check where a line segment intersects the mesh and if it intersects at all
  //
  //! \param[in,out] r0 In: starting position, out: intersection point
  //! \param[in] r1 Ending position
  //! \param[out] ijk Indices of the mesh bin containing the intersection point
  //! \return Whether the line segment connecting r0 and r1 intersects mesh
  bool intersects(Position& r0, Position r1, int* ijk) const;

  // Data members
  xt::xtensor<int, 1> shape_; //!< Number of mesh elements in each dimension

private:
  std::vector<std::vector<double>> grid_;
};

class UnstructuredMeshBase : public Mesh {
public:
  UnstructuredMeshBase(pugi::xml_node node);
  std::string filename_;
};

#ifdef DAGMC

class UnstructuredMesh : public Mesh {

  typedef std::vector<std::pair<double, moab::EntityHandle>> UnstructuredMeshHits;

public:
  UnstructuredMesh() { };
  UnstructuredMesh(pugi::xml_node node);

  //! Determine which bins were crossed by a particle.
  //
  //! \param[in] p Particle to check
  //! \param[out] bins Bins that were crossed
  //! \param[out] lengths Fraction of tracklength in each bin
  void bins_crossed(const Particle* p, std::vector<int>& bins,
                    std::vector<double>& lengths) const;


  bool intersects(Position& r0, Position r1, int* ijk);


private:

//! Finds all intersections with faces of the mesh.
//
//! \param[in] start Staring location
//! \param[in] dir Normalized particle direction
//! \param[in] length of particle track
//! \param[out] Mesh intersections
void
intersect_track(const moab::CartVect& start,
                const moab::CartVect& dir,
                double track_len,
                UnstructuredMeshHits& hits) const;

  //! Calculates the volume for a given tetrahedron handle.
  //
  // \param[in] tet MOAB EntityHandle of the tetrahedron
  double tet_volume(moab::EntityHandle tet) const;

  //! Find the tetrahedron for the given location if
  //! one exists
  //
  //! \param[in]
  //! \return MOAB EntityHandle of tet
  moab::EntityHandle get_tet(const Position& r) const;

  //! Version of get_tet taking Position.
  inline moab::EntityHandle get_tet(const moab::CartVect& r) const {
      return get_tet(Position(r[0], r[1], r[2]));
  };

  //! Check for point containment within a tet, uses
  //! pre-computed barycentric data.
  //
  //! \param[in] r Position to check
  //! \param[in] MOAB terahedron to check
  //! \return True if r is inside, False if r is outside
  bool point_in_tet(const moab::CartVect& r, moab::EntityHandle tet) const;

  //! Compute barycentric coordinate data for all tetrahedra
  //! in the mesh.
  //
  //! \param[in] all_tets MOAB Range of tetrahedral elements
  void compute_barycentric_data(const moab::Range& all_tets);

  //! Translates a MOAB EntityHandle its corresponding bin.
  //
  //! \param[in] eh MOAB EntityHandle to translate
  //! \return Mesh bin
  int get_bin_from_ent_handle(moab::EntityHandle eh) const;

  //! Translates a bin to its corresponding MOAB EntityHandle
  //! for the tetrahedron representing that bin.
  //
  //! \param[in] bin Bin value to translate
  //! \return MOAB EntityHandle of tet
  moab::EntityHandle get_ent_handle_from_bin(int bin) const;

  int get_bin_from_indices(const int* ijk) const override;

  void get_indices(Position r, int* ijk, bool* in_mesh) const override;

  void get_indices_from_bin(int bin, int* ijk) const override;

  std::pair<std::vector<double>, std::vector<double>>
  plot(Position plot_ll, Position plot_ur) const override;

  //! Builds a KDTree for all tetrahedra in the mesh. All
  //! triangles representing 2D faces of the mesh are
  //! added to the tree as well.
  //
  //! \param[in] all_tets MOAB Range of tetrahedra for the tree
  void build_kdtree(const moab::Range& all_tets);

public:
  //! Determine which surface bins were crossed by a particle.
  //
  //! \param[in] p Particle to check
  //! \param[out] bins Surface bins that were crossed
  void surface_bins_crossed(const Particle* p, std::vector<int>& bins) const;

  //! Write mesh data to an HDF5 group.
  //
  //! \param[in] group HDF5 group
  void to_hdf5(hid_t group) const;

  //! Get bin at a given position.
  //
  //! \param[in] r Position to get bin for
  //! \return Mesh bin
  int get_bin(Position r) const;

  std::string get_label_for_bin(int bin) const;

  int n_bins() const override;

  int n_surface_bins() const override;

  std::string filename_; //<! Path to unstructured mesh file

private:
  moab::Range ehs_; //!< Range of tetrahedra EntityHandles in the mesh
  moab::EntityHandle meshset_; //!< Meshset containing all Tets/Tris
  moab::EntityHandle kdtree_root_; //!< Root of the MOAB KDTree
  std::shared_ptr<moab::Interface> mbi_; //!< MOAB instance
  std::unique_ptr<moab::AdaptiveKDTree> kdtree_; //!< MOAB KDTree instance
  std::vector<moab::Matrix3> baryc_data_; //!< Barycentric data for tetrahedra
};

#endif

#ifdef LIBMESH
class LibMesh : public UnstructuredMeshBase {

  typedef std::vector<std::pair<double, const libMesh::Elem*>> UnstructuredMeshHits;

public:
  // constructor
  LibMesh(pugi::xml_node node);

  void bins_crossed(const Particle* p,
                    std::vector<int>& bins,
                    std::vector<double>& lengths) const;

  void intersect_track(libMesh::Point start,
                       libMesh::Point dir,
                       double track_len,
                       UnstructuredMeshHits& hits) const;

  bool elements_share_face(const libMesh::Elem* from,
                           const libMesh::Elem* to,
                           unsigned int side) const;

  int
  get_bin_from_mesh_type(const libMesh::Elem* elem) const;

  int get_bin(Position r) const;

  int n_bins() const;

  void get_indices(Position r, int* ijk, bool* in_mesh) const;

  std::pair<double, const libMesh::Elem*>
  locate_boundary_element(const Position& r0,
                          const Position& r1) const;

  std::pair<double, const libMesh::Elem*>
  locate_boundary_element(const libMesh::Point& start,
                          const libMesh::Point& end) const;


  // triangle intersection methods
  bool plucker_test(std::unique_ptr<const libMesh::Elem> tri,
                    const libMesh::Point& start,
                    const libMesh::Point& dir,
                    double& dist) const;

  double plucker_edge_test(const libMesh::Node& vertexa,
                           const libMesh::Node& vertexb,
                           const libMesh::Point& ray,
                           const libMesh::Point& ray_normal) const;

  double first(const libMesh::Node& a,
               const libMesh::Node& b) const;

  int get_bin_from_indices(const int* ijk) const;

  void get_indices_from_bin(int bin, int* ijk) const;

  libMesh::Elem* get_element_from_bin(int bin) const;

  bool intersects(Position& r0, Position r1, int* ijk) const;

  int n_surface_bins() const;

  void surface_bins_crossed(const Particle* p,
                             std::vector<int>& bins) const;

  std::pair<std::vector<double>, std::vector<double>> plot(Position plot_ll,
                                                           Position plot_ur) const;

  void add_variable(const std::string& var_name);

  void set_variable(const std::string& var_name, int bin, double value);

  void write() const;

  void to_hdf5(hid_t group) const;

private:
  std::unique_ptr<libMesh::Mesh> m_;
  std::unique_ptr<libMesh::PointLocatorBase> point_locator_;
  std::unique_ptr<libMesh::EquationSystems> equation_systems_;
  std::unique_ptr<libMesh::NumericVector<libMesh::Real>> sol_;
  std::map<std::string, unsigned int> variable_map_;
  std::string eq_system_name_;
  libMesh::Elem* first_element_;
  std::set<libMesh::Elem*> boundary_elements_;
};
#endif

//==============================================================================
// Non-member functions
//==============================================================================

//! Read meshes from either settings/tallies
//
//! \param[in] root XML node
void read_meshes(pugi::xml_node root);

//! Write mesh data to an HDF5 group
//
//! \param[in] group HDF5 group
void meshes_to_hdf5(hid_t group);

RegularMesh* get_regular_mesh(int32_t index);

void free_memory_mesh();

} // namespace openmc

#endif // OPENMC_MESH_H
