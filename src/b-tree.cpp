// clang++ -O3 -std=c++17 -march=native standalone.cc
// GCC also compiles, but the performance is slightly worse
// Requires an x86 CPU with at least AVX2

// On Linux, make sure madvise is enabled to use hugepages
// (https://en.algorithmica.org/hpc/cpu-cache/paging/#changing-page-size)

#pragma GCC optimize("O3")
#pragma GCC target("avx2,bmi")

#include <bits/stdc++.h>
#include <x86intrin.h>
#include <sys/mman.h>

// Windows:
// #include "memoryapi.h"

#define _CMP_GT_OQ    0x1e

typedef __m256i reg;
typedef __m256d dreg;

const int N = (1<<16), Q = (1<<22); // <- change these

const int B = 16;
const double INF = std::numeric_limits<double>::max();

constexpr int blocks(int n) {
    return (n + B - 1) / B;
}

constexpr int prev_keys(int n) {
    return (blocks(n) + B) / (B + 1) * B;
}

constexpr int height(int n) {
    return (n <= B ? 1 : height(prev_keys(n)) + 1);
}

constexpr int offset(int h) {
    int k = 0, n = N;
    while (h--) {
        k += blocks(n) * B;
        n = prev_keys(n);
    }
    return k;
}

const int H = height(N), S = offset(H);

double *btree;

void permute(double *node) {
    const reg perm_mask = _mm256_set_epi32(3, 2, 1, 0, 7, 6, 5, 4);
    double* middle = node + 4;
    dreg x = _mm256_loadu_pd(middle);
    x = _mm256_permutevar_pd(x, perm_mask);
    _mm256_storeu_pd(middle, x);
}

void prepare(double *a) {
    const int P = 1 << 21, T = (4 * S + P - 1) / P * P;
    btree = (double*) aligned_alloc(P, T);
    #ifdef __linux__
    madvise(btree, T, MADV_HUGEPAGE);
    #endif

    for (int i = N; i < S; i++)
        btree[i] = INF;

    memcpy(btree, a, 4 * N);

    for (int h = 1; h < H; h++) {
        for (int i = 0; i < offset(h + 1) - offset(h); i++) {
            int k = i / B,
                j = i - k * B;
            k = k * (B + 1) + j + 1;
            for (int l = 0; l < h - 1; l++)
                k *= (B + 1);
            btree[offset(h) + i] = (k * B < N ? btree[k * B] : INF);
        }
    }

    for (int i = offset(1); i < S; i += B)
        permute(btree + i);
}

unsigned direct_rank(dreg x, double* y) {
    dreg a = _mm256_load_pd(y);
    dreg b = _mm256_load_pd(y + 8);

    dreg ca = _mm256_cmp_pd(a, x, _CMP_GT_OQ);
    dreg cb = _mm256_cmp_pd(b, x, _CMP_GT_OQ);

    int mb = _mm256_movemask_pd(cb);
    int ma = _mm256_movemask_pd(ca);

    unsigned mask = (1 << 16);
    mask |= mb << 8;
    mask |= ma;

    return __tzcnt_u32(mask);
}

unsigned permuted_rank(dreg x, double* y) {
    dreg a = _mm256_load_pd(y);
    dreg b = _mm256_load_pd(y + 8);

    dreg ca = _mm256_cmp_pd(a, x, _CMP_GT_OQ);
    dreg cb = _mm256_cmp_pd(b, x, _CMP_GT_OQ);

    dreg c = (dreg) _mm256_packs_epi32((reg)ca, (reg)cb);
    unsigned mask = _mm256_movemask_pd(c);

    return __tzcnt_u32(mask);
}

int lower_bound(double _x) {
    unsigned k = 0;
    dreg x = _mm256_set1_pd(_x - 1);
    for (int h = H - 1; h > 0; h--) {
        unsigned i = permuted_rank(x, btree + offset(h) + k);
        k = k * (B + 1) + (i << 3);
    }
    unsigned i = direct_rank(x, btree + k);
    return btree[k + i];
}

int a[N], q[Q];

// int baseline(int x) {
//     return *std::lower_bound(a, a + N, x);
// }

// double timeit(int (*f)(int)) {
//     clock_t start = clock();

//     int checksum = 0;

//     for (int i = 0; i < Q; i++)
//         checksum ^= f(q[i]);

//     double seconds = double(clock() - start) / CLOCKS_PER_SEC;
//     printf("Checksum: %d\n", checksum);

//     return 1e9 * seconds / Q;
// }

// int main() {
//     printf("N = %d, Q = %d\n", N, Q);

//     std::mt19937 rng(0);

//     for (int i = 0; i < N; i++)
//         a[i] = rng() % (1 << 30);
//     for (int i = 0; i < Q; i++)
//         q[i] = rng() % (1 << 30);

//     a[0] = INF;
//     std::sort(a, a + N);
//     prepare(a);

//     double x = timeit(baseline);
//     double y = timeit(lower_bound);

//     printf("std::lower_bound: %.2f\n", x);
//     printf("S+ tree: %.2f\n", y);

//     printf("Speedup: %.2f\n", x / y);

//     return 0;
// }