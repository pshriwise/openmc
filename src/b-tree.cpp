#include <iostream>

#include <x86intrin.h>
#include <limits>

#pragma GCC optimize("O3")
#pragma GCC target("avx2,bmi")

#include <bits/stdc++.h>
#include <x86intrin.h>
#include <sys/mman.h>

typedef __m256i reg;
typedef __m256d dreg;

const int B = 8;
const double INF = std::numeric_limits<double>::max();

// state variables
int N = (1<<20);

int nblocks;
double *_a;
double *btree;

int blocks(int n) {
    return (n + B - 1) / B;
}

int prev_keys(int n) {
    return (blocks(n) + B) / (B + 1) * B;
}

int height(int n) {
    return (n <= B ? 1 : height(prev_keys(n)) + 1);
}

int offset(int h) {
    int k = 0, n = N;
    while (h--) {
        k += blocks(n) * B;
        n = prev_keys(n);
    }
    return k;
}

int H = height(N), S = offset(H);

void prepare(double *a, int _n) {
    const int P = 1 << 21, T = (8 * S + P - 1) / P * P;
    btree = (double*) std::aligned_alloc(P, T);
    madvise(btree, T, MADV_HUGEPAGE);

    for (int i = N; i < S; i++) btree[i] = INF;

    memcpy(btree, a, 8*N);

    // create internal nodes
    for (int h = 1; h < H; h++) {
        for (int i = 0; i < offset(h + 1) - offset(h); i++) {
            int k = i / B,
                j = i - k * B;
            k = k * (B + 1) + j + 1; // compare right
            // and then always to the left
            for (int l = 0; l < h - 1; l++)
                k *= (B + 1);
            btree[offset(h) + i] = (k * B < N ? btree[k * B] : INF);
        }
    }
}

int direct_rank(dreg x, double* y_ptr) {
    dreg y1 = _mm256_load_pd(y_ptr);
    dreg y2 = _mm256_load_pd(y_ptr + 4);

    dreg mask1 = _mm256_cmp_pd(x, y1, _CMP_GT_OQ);
    dreg mask2 = _mm256_cmp_pd(x, y2, _CMP_GT_OQ);

    int mask = ~(_mm256_movemask_pd(mask1) + (_mm256_movemask_pd(mask2) << 4));

    return __builtin_ffs(mask) - 1;
}

double lower_bound(double _x) {
    unsigned k = 0;
    dreg x_vec = _mm256_set1_pd(_x);
    for (int h = H - 1; h > 0; h--) {
       int i = direct_rank(x_vec, btree + offset(h) + k);
       k = k * (B + 1) + (i << 3);
    }
    int i = direct_rank(x_vec, btree + k);
    return btree[k + i];
}

#ifdef BTREE_BENCHMARK

int main(int argc, char* argv[]) {
    int n = (argc > 1 ? atoi(argv[1]) : N);
    int m = (argc > 2 ? atoi(argv[2]) : 1<<20);

    N = n;
    H = height(N);
    S = offset(H);

    double *a = new double[n];
    double *q = new double[m];

    // std::cout << "Array size: " << n << std::endl;
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

    prepare(a, n);

    int checksum = 0;
    clock_t start = clock();

    for (int i = 0; i < m; i++) {
        // std::cout << "Query: " << q[i] << std::endl;
        // double result = lower_bound(q[i]);
        // std::cout << "Result: " << result << std::endl;
        checksum ^= (int)lower_bound(q[i]);
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
        checksum ^= (int)*std::lower_bound(a, a+n, q[i]);
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

