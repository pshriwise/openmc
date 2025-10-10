#ifndef OPENMC_TALLIES_FILTER_UNIVERSE_CELL_H
#define OPENMC_TALLIES_FILTER_UNIVERSE_CELL_H

#include <cstdint>
#include <unordered_map>

#include <gsl/gsl-lite.hpp>

#include "openmc/tallies/filter.h"
#include "openmc/tallies/filter_cell.h"
#include "openmc/vector.h"

namespace openmc {

//==============================================================================
//! Specifies which geometric universes tally events reside in.
//==============================================================================

class UniverseCellFilter : public CellFilter {
public:
  //----------------------------------------------------------------------------
  // Constructors, destructors

  ~UniverseCellFilter() = default;

  //----------------------------------------------------------------------------
  // Methods
  void from_xml(pugi::xml_node node) override;

  //----------------------------------------------------------------------------
  // Accessors
};

} // namespace openmc
#endif // OPENMC_TALLIES_FILTER_UNIVERSE_CELL_H
