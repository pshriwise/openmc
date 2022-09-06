
#ifndef OPENMC_INTERPOLATE_H
#define OPENMC_INTERPOLATE_H

#include <cmath>
#include <type_traits>
#include <vector>

#include "openmc/error.h"
#include "openmc/search.h"

namespace openmc {

template<class itx>
struct Interpolator {

  Interpolator(itx arr_begin, itx arr_end, double x, size_t idx,
    Interpolation i = Interpolation::lin_lin)
    : interpolation_(i), idx_(idx)
  {
    set_factor(arr_begin + idx_, x);
  }

  Interpolator(itx arr_begin, itx arr_end, double x,
    Interpolation i = Interpolation::lin_lin)
    : interpolation_(i)
  {
    // set index into array
    if (x < *arr_begin) {
      idx_ = 0;
      interpolation_factor_ = 0.0;
    } else if (x > *(arr_end - 1)) {
      idx_ = arr_end - arr_begin - 2;
      interpolation_factor_ = 1.0;
    } else {
      idx_ = lower_bound_index(arr_begin, arr_end, x);
      set_factor(arr_begin + idx_, x);
    }
  }

  void set_factor(itx pos, double x)
  {
    double x0 = *pos;
    double x1 = *(pos + 1);

    // compute interpolation factor
    switch (interpolation_) {
    case Interpolation::lin_lin:
    case Interpolation::log_lin:
      interpolation_factor_ = (x - x0) / (x1 - x0);
      break;
    case Interpolation::lin_log:
    case Interpolation::log_log:
      interpolation_factor_ = log(x / x0) / log(x1 / x0);
      break;
    default:
      fatal_error("Unrecognized interpolation");
      break;
    }
  }

  // performs interpolation for on specified y values
  double operator()(double y0, double y1)
  {
    switch (interpolation_) {
    case Interpolation::lin_lin:
    case Interpolation::lin_log:
      return y0 + interpolation_factor_ * (y1 - y0);
    case Interpolation::log_lin:
    case Interpolation::log_log:
      return y0 * exp(interpolation_factor_ * log(y1 / y0));
    default:
      fatal_error("Unrecognized interpolation");
    }
  }

  template<class ity>
  double operator()(ity arr_begin)
  {
    double y0 = *(arr_begin + idx_);
    double y1 = *(arr_begin + idx_ + 1);

    return (*this)(y0, y1);
  }

  double operator()(const std::vector<double>& ys)
  {
    return (*this)(ys.begin());
  }

  // Accessors
  size_t idx() const { return idx_; }
  double factor() const { return interpolation_factor_; }

  // Data members
  Interpolation interpolation_ {0.0};
  itx begin_;
  itx end_;
  size_t idx_ {0};
  double interpolation_factor_;
};

inline double interpolate_lagrangian(const std::vector<double>& xs,
  const std::vector<double>& ys, double x, int order)
{
  int idx = lower_bound_index(xs.begin(), xs.end(), x);

  std::vector<double> coeffs;

  for (int i = 0; i < order + 1; i++) {
    double numerator {1.0};
    double denominator {1.0};
    for (int j = 0; j < order; j++) {
      if (i == j)
        continue;
      numerator *= (x - xs[idx + j]);
      denominator *= (xs[idx + i] - xs[idx + j]);
    }
    coeffs.push_back(numerator / denominator);
  }

  return std::inner_product(
    coeffs.begin(), coeffs.end(), ys.begin() + idx, 0.0);
}

// Pseudo-constructor for Interpolator class to handle CTAD, constructor can be
// called directly when we move to C++17
template <typename F>
Interpolator<F> FixedInterpolator(F&& begin, F&& end, double x, Interpolation i = Interpolation::lin_lin) {
    return Interpolator<F>{std::forward<F>(begin), std::forward<F>(end), x, i};
}

} // namespace openmc

#endif