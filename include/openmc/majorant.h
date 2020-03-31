//! \file majorant.h
//! \brief Majorant cross section type

#ifndef OPENMC_MAJORANT_H
#define OPENMC_MAJORANT_H

#include<gsl/gsl>
#include <vector>

#include "openmc/settings.h"

namespace openmc {

class Majorant {

  struct EnergyGrid {
    std::vector<int> grid_index;
    std::vector<double> energy;
  };

  struct XS {
    XS(std::vector<double> energies,
       std::vector<double> total_xs)
    : energies_(energies), total_xs_(total_xs), idx_(0)
    { Expects(energies_.size() == total_xs_.size()); }

    //! \brief Return the current energy and total cross section values
    std::pair<double, double> get() const;

    //! \brief Return the energy of the current index
    double get_e() const;

    //! \brief Return the cross section of the current index
    double get_xs() const;

    //! \brief Return the current energy value and advance one
    double pop_e();

    //! \brief Advance the index until past the input energy value
    void advance(double energy);

    //! \brief Indicate if all values in the cross section have been visited
    bool complete() const;

    //! \brief Return the previous energy, cross section value pair
    std::pair<double, double> prev() const;

    //! \brief Return the previous energy value
    double prev_e() const;

    //! \brief Return the previous cross section value
    double prev_xs() const;

    //! \brief Increment the cross section index by i
    inline
    XS& operator +=(size_t i) { this->idx_ += i; return *this; }

    //! \brief Increment the cross section index by one
    inline
    XS& operator ++(int i) { this->idx_ += 1; return *this; }

    std::vector<double> energies_;
    std::vector<double> total_xs_;
    size_t idx_;
  };

 public:
  void write_ascii() const;

  //! \brief Initialize the energy grid mapping for the majorant xs
  void init_grid();

  //! \brief Update the majorant using values from another cross section
  void update(std::vector<double> energies_other,
              std::vector<double> xs_other);

    // data members
 public:
  std::vector<int> nuclides; // index of nuclides applied
  std::vector<double> xs_; // cross section values
  EnergyGrid grid_;
  constexpr static double safety_factor {1.01};
}; // class Majorant

  bool intersect_2D(std::pair<double, double> p1,
                    std::pair<double, double> p2,
                    std::pair<double, double> p3,
                    std::pair<double, double> p4,
                    std::pair<double, double>& out);

  bool is_above(std::pair<double, double> p1,
                std::pair<double, double> p2,
                std::pair<double, double> p3);

  void create_majorant();
}

#endif // OPENMC_MAJORANT_H
