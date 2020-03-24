
#include <ostream>

#include "openmc/majorant.h"
#include "openmc/nuclide.h"

namespace openmc {

  void Majorant::write_ascii() const {

    std::ofstream of("majorant.txt");

    for (int i = 0; i < xs_.size(); i++) {
      of << e_[i] << "\t" << xs_[i];
    }

    of.close();
  }

}
