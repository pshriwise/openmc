#include "openmc/tallies/filter_universe_cell.h"

#include <fmt/core.h>

#include "openmc/cell.h"
#include "openmc/error.h"
#include "openmc/xml_interface.h"

namespace openmc {

void UniverseCellFilter::from_xml(pugi::xml_node node)
{
  // Get material IDs and convert to indices in the global materials vector
  auto universes = get_node_array<int32_t>(node, "bins");
  for (auto& u : universes) {
    auto search = model::universe_map.find(u);
    if (search == model::universe_map.end()) {
      throw std::runtime_error {fmt::format(
        "Could not find universe {} specified on tally filter.", u)};
    }
    u = search->second;
  }

  vector<int32_t> cell_ids;

  for (auto u : universes) {
    for (auto c : model::universes[u]->cells_) {
      cell_ids.push_back(c);
    }
  }

  this->set_cells(cell_ids);
}

} // namespace openmc
