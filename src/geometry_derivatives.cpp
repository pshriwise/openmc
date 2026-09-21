#include "openmc/geometry_derivatives.h"
#include "openmc/error.h"

namespace openmc {

namespace model {
  std::unordered_map<int32_t, int32_t> geometry_derivatives_map;
  vector<unique_ptr<GeometryDerivative>> geometry_derivatives;
}

GeometryDerivative::GeometryDerivative(pugi::xml_node node)
{
  if (!check_for_node(node, "id")) {
    fatal_error("Must specify id of geometry derivative in geometry_derivatives XML file.");
  }
  id_ = std::stoi(get_node_value(node, "id"));

  if (!check_for_node(node, "tally_id")) {
    fatal_error("Must specify tally_id of geometry derivative in geometry_derivatives XML file.");
  }
  tally_id_ = std::stoi(get_node_value(node, "tally_id"));

  if (!check_for_node(node, "cell_id")) {
    fatal_error("Must specify cell_id of geometry derivative in geometry_derivatives XML file.");
  }
  cell_id_ = std::stoi(get_node_value(node, "cell_id"));
}

void read_geometry_derivatives(pugi::xml_node node)
{
  if (!check_for_node(node, "geometry_derivatives")) {
    return;
  }

  for (pugi::xml_node geom_deriv_node : node.children("geometry_derivative")) {
    model::geometry_derivatives.push_back(std::make_unique<GeometryDerivative>(geom_deriv_node));
    model::geometry_derivatives_map[model::geometry_derivatives.back()->id()] = model::geometry_derivatives.size() - 1;
  }
}

} // namespace openmc