
#ifndef OPENMC_INTERPOLATE_H
#define OPENMC_INTERPOLATE_H

#include <cmath>
#include <type_traits>
#include <vector>

#include "openmc/error.h"
#include "openmc/search.h"

namespace openmc {

struct FixedInterpolator {
  FixedInterpolator(const std::vector<double>& xs, double x,
    Interpolation i = Interpolation::lin_lin)
    : interpolation_(i)
  {
    // set index into array
    if (x < xs.front())
      idx_ = 0;
    else if (x > xs.back())
      idx_ = xs.size() - 2;
    else
      idx_ = lower_bound_index(xs.begin(), xs.end(), x);

    set_factor(xs.begin() + idx_, x);
  }

  FixedInterpolator(const std::vector<double>& xs, double x, size_t idx,
    Interpolation i = Interpolation::lin_lin)
    : interpolation_(i), idx_(idx)
  {
    set_factor(xs.begin() + idx_, x);
  }

  template<class It>
  FixedInterpolator(It arr_begin, It arr_end, double x, size_t idx,
    Interpolation i = Interpolation::lin_lin)
    : interpolation_(i), idx_(idx)
  {
    set_factor(arr_begin + idx_, x);
  }

  template<class It>
  FixedInterpolator(It arr_begin, It arr_end, double x,
    Interpolation i = Interpolation::lin_lin)
    : interpolation_(i)
  {
    // set index into array
    if (x < *arr_begin)
      idx_ = 0;
    else if (x > *(arr_end - 1))
      idx_ = arr_end - arr_begin - 2;
    else
      idx_ = lower_bound_index(arr_begin, arr_end, x);

    set_factor(arr_begin + idx_, x);
  }

  template<class It>
  void set_factor(It pos, double x)
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

  template<class It>
  double operator()(It arr_begin)
  {
    double y0 = *(arr_begin + idx_);
    double y1 = *(arr_begin + idx_ + 1);

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

  double operator()(const std::vector<double>& ys)
  {
    return (*this)(ys.begin());
  }

  // Accessors
  size_t idx() const { return idx_; }
  double factor() const { return interpolation_factor_; }

  // Data members
  Interpolation interpolation_;
  size_t idx_;
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

} // namespace openmc

#endif