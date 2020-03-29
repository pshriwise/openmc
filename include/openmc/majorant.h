//! \file majorant.h
//! \brief Majorant cross section type

#ifndef OPENMC_MAJORANT_H
#define OPENMC_MAJORANT_H

#include<gsl/gsl>
#include <vector>

namespace openmc {

class Majorant {

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

    //! \brief Advance the index until past the input energy value
    void advance(double energy);

    //! \brief Indicate if all values in the cross section have been visited
    bool complete() const;

    //! \brief Increment the cross section index
    // inline
    // XS& operator +=(size_t i) { this->idx_ += i; return *this; }

    std::vector<double> energies_;
    std::vector<double> total_xs_;
    size_t idx_;
  };

 public:
  void write_ascii() const;

    // data members
 public:
  std::vector<int> nuclides; // index of nuclides applied
  std::vector<double> e_; // energy points
  std::vector<double> xs_; // cross section values
}; // class Majorant

  static std::pair<double, double>
    intersect_2D(double x1, double y1,
                 double x2, double y2,
                 double x3, double y3,
                 double x4, double y4);



}

#endif // OPENMC_MAJORANT_H
