
#include <bits/stdc++.h>
#include <cstdlib>
#include <limits>
#include <vector>
#include <sys/mman.h>
#include <x86intrin.h>

#include "openmc/constants.h"

#pragma GCC optimize("O3")
#pragma GCC target("avx2,bmi")

namespace openmc {

typedef __m256i reg;
typedef __m256d dreg;

class SearchArray {

  // static b+ tree branching ratio
  static const int B_ = 8;

public:
  SearchArray(const std::vector<double>& arr);

  SearchArray(const double* arr, int N);

  ~SearchArray() { free(btree); }

private:

  inline int blocks(int n) const {
    return (n + B_ - 1) / B_;
  }

  inline int prev_keys(int n) const {
    return (blocks(n) + B_) / (B_ + 1) * B_;
  }

  inline int height(int n) const {
    return (n <= B_ ? 1 : height(prev_keys(n)) + 1);
  }

  int offset(int h) const {
    int k = 0, n = N_;
    while (h--) {
        k += blocks(n) * B_;
        n = prev_keys(n);
    }
    return k;
  }

  static int direct_rank(dreg x, double* y_ptr) {
    dreg y1 = _mm256_load_pd(y_ptr);
    dreg y2 = _mm256_load_pd(y_ptr + 4);

    dreg mask1 = _mm256_cmp_pd(x, y1, _CMP_GT_OQ);
    dreg mask2 = _mm256_cmp_pd(x, y2, _CMP_GT_OQ);

    int mask = ~(_mm256_movemask_pd(mask1) + (_mm256_movemask_pd(mask2) << 4));

    return __builtin_ffs(mask) - 1;
  }

  void prepare(const double* arr);

public:
  double lower_bound(double _x) const;

  // data members
  private:
      const int N_;
      const int H_;
      double* btree;
};

} // namespace openmc