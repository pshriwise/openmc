
#include <fstream>

#include "openmc/majorant.h"
#include "openmc/nuclide.h"

namespace openmc {

  static std::pair<double, double>
  intersect_2D(double x1, double y1,
               double x2, double y2,
               double x3, double y3,
               double x4, double y4) {

    double denominator = (x4 - x3) * (y1 - y2) - (x1 - x2) * (y4 - y3);

    // if the lines are parallel, return no intersection
    if (denominator == 0.0) { return {-1.0, -1.0}; }

    double numerator = (y3 - y4) * (x1 - x3) + (x4 - x3) * (y1 - y3);

    double t = numerator / denominator;

    if (0.0 <= t and t <= 1.0) {
      double x = x1 + (x2 - x1) * t;
      double m = (y2 - y1) / (x2 - x1);
      double y = y1 + (x - x1) * m;

      return {x, y};
    }

    return {-1, -1};
  }


  // Majorant method definitions
  void Majorant::write_ascii() const {

    std::ofstream of("majorant.txt");

    for (int i = 0; i < xs_.size(); i++) {
      of << e_[i] << "\t" << xs_[i];
    }

    of.close();
  }

  // Majorant XS definitions

  std::pair<double, double>
  Majorant::XS::get() const { return {energies_[idx_], total_xs_[idx_]}; }

  double Majorant::XS::get_e() const { return energies_[idx_]; }

  double Majorant::XS::get_xs() const { return total_xs_[idx_]; }

  void Majorant::XS::advance(double energy) {
    double e = energies_[idx_];
    while (e <= energy && !this->complete()) { idx_++; }
  }

  bool Majorant::XS::complete() const { return idx_ >= energies_.size(); }

}
