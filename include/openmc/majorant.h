//! \file majorant.h
//! \brief Majorant cross section type

#ifdef OPENMC_PARTICLE_H
#define OPENMC_PARTICLE_H

#include <vector>

namespace openmc {

class Majorant {

 public:
  void write_ascii() const;

    // data members
 public:
  std::vector<int> nuclides; // index of nuclides applied
  std::vector<double> e_; // energy points
  std::vector<double> xs_; // cross section values
};

}

#endif // OPENMC_PARTICLE_H
