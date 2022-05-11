#include <iostream>

#include "openmc/b-tree.h"

namespace openmc {

SearchArray::SearchArray(const std::vector<double>& arr) : N_(arr.size()), H_(height(N_)) {
  prepare(arr.data());
};

SearchArray::SearchArray(const double* arr, int N) : N_(N), H_(height(N_)) {
  prepare(arr);
};

double SearchArray::lower_bound(double _x) const {
  unsigned k = 0;
  dreg x_vec = _mm256_set1_pd(_x);
  for (int h = H_ - 1; h > 0; h--) {
    int i = direct_rank(x_vec, btree + offset(h) + k);
    k = k * (B_ + 1) + (i << 3);
  }
  int i = direct_rank(x_vec, btree + k);
  return btree[k + i];
}

void SearchArray::prepare(const double* arr) {
  int S = offset(H_);
  // compute size of array in bytes
  const int P = 1 << 21;
  const int T = (8 * S + P - 1) / P * P;
  btree = (double*) aligned_alloc(P, T);
  madvise(btree, T, MADV_HUGEPAGE);

  for (int i = N_; i < S; i++) btree[i] = INFTY;
  memcpy(btree, arr, 8*N_);

  // create internal nodes
  for (int h = 1; h < H_; h++) {
    for (int i = 0; i < offset(h + 1) - offset(h); i++) {
      int k = i / B_,
        j = i - k * B_;
      k = k * (B_ + 1) + j + 1; // compare right
      // and then always to the left
      for (int l = 0; l < h - 1; l++)
        k *= (B_ + 1);
      btree[offset(h) + i] = (k * B_ < N_ ? btree[k * B_] : INFTY);
    }
  }
}

}

#ifdef BTREE_BENCHMARK

int main(int argc, char* argv[]) {
    int n = (argc > 1 ? atoi(argv[1]) : 1<<20);
    int m = (argc > 2 ? atoi(argv[2]) : 1<<21);

    // N = n;
    // H = height(N);
    // S = offset(H);


    double *a = new double[n];
    double *q = new double[m];

    std::cout << "Array size: " << n << std::endl;
    for (int i = 0; i < n; i++)
        a[i] = (double)rand() / (double)RAND_MAX;

    // std::cout << "# Queries: " << m << std::endl;
    for (int i = 0; i < m; i++) {
        q[i] = (double)rand() / (double)RAND_MAX;
    }

    // std::cout << "Tree height: " << H << std::endl;

    std::cout << std::endl;

    a[0] = RAND_MAX;
    std::sort(a, a + n);

    openmc::SearchArray arr{a, n};

    int checksum = 0;
    clock_t start = clock();

    for (int i = 0; i < m; i++) {
        // std::cout << "Query: " << q[i] << std::endl;
        // double result = lower_bound(q[i]);
        double result = arr.lower_bound(q[i]);
        checksum ^= (int)result;

        // std::cout << "Result: " << result << std::endl;
    }

    float b_seconds =  float(clock() - start) / CLOCKS_PER_SEC;


    // std::cout << "B+ Tree" << std::endl;
    // std::cout << "--------" << std::endl;
    // printf("%.2f ns per query\n", 1e9 * b_seconds / m);
    std::cerr << "Checksum: " << checksum << std::endl;
    // std::cout << "--------" << std::endl;
    // std::cout << std::endl;

    checksum = 0;
    start = clock();

    for (int i = 0; i < m; i++) {
      double result = *std::lower_bound(a, a+n, q[i]);
      checksum ^= (int)result;
    }

    float std_seconds = float(clock() - start) / CLOCKS_PER_SEC;

    // std::cout << "std::lower_bound" << std::endl;
    // std::cout << "-------" << std::endl;
    // printf("%.2f ns per query\n", 1e9 * std_seconds / m);
    // std::cout << "Checksum: " << checksum << std::endl;
    // std::cout << "-------" << std::endl;

    // std::cout << std::endl;
    // std::cout << "Speedup: " << std_seconds / b_seconds << std::endl;

    std::cout << n << ", " << 1e9 * b_seconds / m << ", " << 1e9 * std_seconds / m << ", " << checksum << std::endl;

    return 0;
}

#endif

