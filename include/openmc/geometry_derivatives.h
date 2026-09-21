#include <unordered_map>

#include "openmc/memory.h"
#include "openmc/vector.h"
#include "openmc/xml_interface.h"

namespace openmc {

class GeometryDerivative;

namespace model {
  extern std::unordered_map<int32_t, int32_t> geometry_derivatives_map;
  extern vector<unique_ptr<GeometryDerivative>> geometry_derivatives;
}

void read_geometry_derivatives(pugi::xml_node node);

class GeometryDerivative {

  public:
  // Constructors
  GeometryDerivative(pugi::xml_node node);

  // Accessors
  int32_t id() const { return id_; }
  int32_t tally_id() const { return tally_id_; }
  int32_t cell_id() const { return cell_id_; }

  private:
  int32_t id_;
  int32_t tally_id_;
  int32_t cell_id_;
};

} // namespace openmc