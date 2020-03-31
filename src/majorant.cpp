
#include <fstream>

#include "openmc/constants.h"
#include "openmc/majorant.h"
#include "openmc/nuclide.h"

#include "xtensor/xview.hpp"

namespace openmc {

  void create_majorant() {
    // create a majorant XS
    Majorant majorant;

    for (const auto& nuclide : data::nuclides) {
      for (int t = 0; t < nuclide->kTs_.size(); t++) {
        auto total = xt::view(nuclide->xs_[t], xt::all(), 0);
        std::vector<double> xs;
        for (auto val : total) { xs.push_back(val); }
        auto energies = nuclide->grid_[t].energy;
        majorant.update(energies, xs);

        // write nuclide to file
        auto T = nuclide->kTs_[t];
        std::ofstream of(nuclide->name_ + "_" + std::to_string(T/K_BOLTZMANN) + ".txt");
        for (int i = 0; i < energies.size(); i++) {
          of << energies[i] << "\t" << xs[i] << "\n";
        }
        of.close();
      }
    }

    majorant.init_grid();
    majorant.write_ascii();
  }

  bool intersect_2D(std::pair<double, double> p1,
                    std::pair<double, double> p2,
                    std::pair<double, double> p3,
                    std::pair<double, double> p4,
                    std::pair<double, double>& out) {

    double denominator = (p4.first - p3.first) * (p1.second - p2.second) -
                         (p1.first - p2.first) * (p4.second - p3.second);

    // if the lines are parallel, return no intersection
    if (denominator == 0.0) { return false; }

    double numerator = (p3.second - p4.second) * (p1.first - p3.first) +
                       (p4.first - p3.first) * (p1.second - p3.second);

    double t = numerator / denominator;

    if (t < 0.0 || t > 1.0) { return false; }

    double x = p1.first + (p2.first - p1.first) * t;
    double m = (p2.second - p1.second) / (p2.first - p1.first);
    double y = p1.second + (x - p1.first) * m;
    out = {x, y};
    return true;
  }

  bool is_above(std::pair<double, double> p1,
                std::pair<double, double> p2,
                std::pair<double, double> p3) {
    // if the line is vertical, use
    // comparison of x values
    if (fabs(p2.first - p1.first) < FP_COINCIDENT) {
      return p3.first < p2.first;
    } else {
      double slope = (p2.second - p1.second) / (p2.first - p1.first);
      double val = p1.second + slope * (p3.first - p1.first);
      return val < p3.second;
    }

  }

  // Majorant method definitions
  void Majorant::write_ascii() const {

    std::ofstream of("majorant.txt");

    for (int i = 0; i < xs_.size(); i++) {
      of << grid_.energy[i] << "\t" << xs_[i] << "\n";
    }

    of.close();
  }

void Majorant::init_grid() {
  int neutron = static_cast<int>(Particle::Type::neutron);
  double E_min = data::energy_min[neutron];
  double E_max = data::energy_max[neutron];
  int M = settings::n_log_bins;

  // Determine equal-logarithmic energy spacing
  double spacing = std::log(E_max/E_min)/M;

  // Create equally log-spaced energy grid
  auto umesh = xt::linspace(0.0, M*spacing, M+1);

  // Resize array for storing grid indices
  grid_.grid_index.resize(M + 1);

  // Determine corresponding indices in nuclide grid to energies on
  // equal-logarithmic grid
  int j = 0;
  for (int k = 0; k <= M; ++k) {
    while (std::log(grid_.energy[j + 1]/E_min) <= umesh(k)) {
      // Ensure that for isotopes where maxval(grid.energy) << E_max that
      // there are no out-of-bounds issues.
      if (j + 2 == grid_.energy.size()) break;
      ++j;
    }
    grid_.grid_index[k] = j;
  }
}

  void Majorant::update(std::vector<double> energy_other,
                        std::vector<double> xs_other) {

    // remove the additional 5 percent from the current xs values
    // std::transform(xs_.begin(), xs_.end(), xs_.begin(),[](double xs_val) { return xs_val / Majorant::safety_factor; });

    XS xs_a(grid_.energy, xs_);
    XS xs_b(energy_other, xs_other);

    // early exit checks
    if (xs_b.complete()) { return; }

    if (xs_a.complete()) {
      grid_.energy = energy_other;
      xs_ = xs_other;
      return;
    }

    // resulting output of this algorithm
    std::vector<double> e_out, xs_out;
    std::vector<size_t> mask;

    // references we'll use to keep track of which
    // is the higher-value cross section
    XS& current_xs = xs_a;
    XS& other_xs = xs_b;

    // if the other cross section starts at a lower energy, it's
    // value is considered to be higher
    if (other_xs.get_e() < current_xs.get_e())
      std::swap(current_xs, other_xs);

    // if the two cross sections start at the same energy
    // pick the one with the higher xs value
    if (other_xs.get_e() == current_xs.get_e() && other_xs.get_xs() > current_xs.get_xs())
      std::swap(current_xs, other_xs);

    // add the first point to the final cross section
    e_out.push_back(current_xs.get_e());
    xs_out.push_back(current_xs.get_xs());
    current_xs++;

    // start looping over values
    while (true) {

      if (current_xs.complete() || other_xs.complete()) { break; }

      auto current_next = current_xs.get();
      auto other_next = other_xs.get();
      std::pair<double, double> last_pnt = {e_out.back(), xs_out.back()};

      // the next point added to the cross section is selected based 4 different cases
      //
      // above: the other cross section point above the line between the last point added and the current xs's next point
      // below: ther other cross section point is below the line betwee the last point added and the current xs's next point
      // nearer: the other cross section point is closer in energy to the last point added
      // farther: the other cross section point is farther in energy from the last point added
      //
      // Case 1: other xs value is above and nearer
      //
      //         Check for an intersection with the line segment between last_pnt and current_next.
      //         If an intersection exists, insert a point in the majorant at that location and
      //         swap the current and other cross sections.
      //
      // Case 2: other xs value is above and farther
      //
      //         Check for an intersection with the line segment between last_pnt and current_next.
      //         Apply the same treatment as Case 1.
      //
      // Case 3: other xs value is below and nearer
      //
      //         Insert a point at the energy of the other xs value on the line segment between
      //         last_pnt and current_next. This point is superfluous and purely for the algorithm
      //         robustness, so the index of this point is added to the mask of indices to be removed
      //         later.
      //
      // Case 4: other xs value is below and farther
      //
      //         Insert the location of current_next to the output cross section and continue

      bool above = is_above(last_pnt, current_next, other_next);
      bool nearer = other_next.first < current_next.first;

      // for clarity
      bool below = !above;
      bool farther = !nearer;

      // Cases 1 & 2
      if (above) {
        // setup intersection check
        auto p1 = last_pnt;
        auto p2 = current_next;
        auto p3 = other_xs.prev();
        auto p4 = other_next;

        std::pair<double, double> intersection;
        auto intersection_found = intersect_2D(p1, p2, p3, p4, intersection);

        if (intersection_found) {
          e_out.push_back(intersection.first);
          xs_out.push_back(intersection.second);
          std::swap(current_xs, other_xs);
        } else {
          e_out.push_back(current_next.first);
          xs_out.push_back(current_next.second);
        }
      // Case 3
      } else if (below and nearer) {
        // compute the y value
        double slope = (current_next.second - last_pnt.second) /
                       (current_next.first - last_pnt.first);
        double xs_val = last_pnt.second + slope * (other_next.first - last_pnt.first);
        e_out.push_back(other_next.first);
        xs_out.push_back(xs_val);
        // add index of this point to the mask
        mask.push_back(e_out.size() - 1);
      // Case 4
      } else if (below and farther) {
        e_out.push_back(current_next.first);
        xs_out.push_back(current_next.second);
      }

      current_xs.advance(e_out.back());
      other_xs.advance(e_out.back());
    }

    // one or both of the cross sections should be complete
    Expects(xs_a.complete() || xs_b.complete());

    // remove any superfluous points
    std::sort(mask.begin(), mask.end(), std::greater<double>());
    for (auto idx : mask) {
      e_out.erase(e_out.begin() + idx);
      xs_out.erase(xs_out.begin() + idx);
    }

    // add any additional data to the output xs
    while (!xs_a.complete()) {
      e_out.push_back(xs_a.get_e());
      xs_out.push_back(xs_a.get_xs());
      xs_a++;
    }

    while (!xs_b.complete()) {
      e_out.push_back(xs_b.get_e());
      xs_out.push_back(xs_b.get_xs());
      xs_b++;
    }

    // increase all xs values in the majorant by 5 percent
    // std::transform(xs_out.begin(), xs_out.end(), xs_out.begin(), [](double xs_val) { return xs_val * safety_factor; });

    // update the values of the majorant
    grid_.energy = std::move(e_out);
    xs_ = std::move(xs_out);

  }

  // Majorant XS definitions

  std::pair<double, double>
  Majorant::XS::get() const { return {energies_.at(idx_), total_xs_.at(idx_)}; }

  double Majorant::XS::get_e() const { return energies_.at(idx_); }

  double Majorant::XS::get_xs() const { return total_xs_.at(idx_); }

  std::pair<double, double>
  Majorant::XS::prev() const { return {energies_.at(idx_ - 1), total_xs_.at(idx_ - 1)}; }

  double Majorant::XS::prev_e() const { return energies_.at(idx_ - 1); }

  double Majorant::XS::prev_xs() const { return total_xs_.at(idx_ - 1); }

  void Majorant::XS::advance(double energy) {
    double e = energies_[idx_];
    while (e <= energy && !this->complete()) { e = energies_[++idx_]; }
  }

  bool Majorant::XS::complete() const { return idx_ >= energies_.size(); }



}
