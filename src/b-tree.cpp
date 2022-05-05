#include <iostream>

#include <x86intrin.h>
#include <limits>

#pragma GCC optimize("O3")
#pragma GCC target("avx2,bmi")

#include <bits/stdc++.h>
#include <x86intrin.h>
#include <sys/mman.h>

#ifdef BTREE_BENCHMARK
#include <bits/stdc++.h>

#ifndef N
#define N (1<<20)
#endif

void prepare(double *a, int n);
double lower_bound(double x);

int main(int argc, char* argv[]) {
    int n = (argc > 1 ? atoi(argv[1]) : N);
    int m = (argc > 2 ? atoi(argv[2]) : 1<<20);

    double *a = new double[n];
    double *q = new double[m];

    std::cout << "Array size: " << n << std::endl;
    for (int i = 0; i < n; i++)
        a[i] = (double)rand() / (double)RAND_MAX;

    std::cout << "Queries: " << m << std::endl;
    for (int i = 0; i < m; i++) {
        q[i] = (double)rand() / (double)RAND_MAX;
        //std::cout << "Query val: " << q[i] << std::endl;
    }


    a[0] = RAND_MAX;
    std::sort(a, a + n);

    prepare(a, n);

    int checksum = 0;
    clock_t start = clock();

    for (int i = 0; i < m; i++) {
        // std::cout << i << std::endl;
        // std::cout << "Query: " << q[i] << std::endl;
        // double result = lower_bound(q[i]);
        // std::cout << "Result: " << result << std::endl;
        checksum ^= (int)lower_bound(q[i]);
        // checksum ^= (int)*std::lower_bound(a, a+n, q[i]);
    }

    float seconds = float(clock() - start) / CLOCKS_PER_SEC;

    printf("%.2f ns per query\n", 1e9 * seconds / m);
    std::cout << "Checksum: " << checksum << std::endl;

    return 0;
}

#endif

typedef __m256i reg;
typedef __m256d dreg;

const int B = 8;
const double INF = std::numeric_limits<double>::max();

int n;
int nblocks;
double *_a;
double (*btree)[B];

int go(int k, int i) { return k * (B + 1) + i + 1; }

void build(int k = 0) {
    static int t = 0;
    if (k < nblocks) {
        for (int i = 0; i < B; i++) {
            build(go(k, i));
            btree[k][i] = (t < n ? _a[t++] : INF);
        }
        build(go(k, B));
    }
}

void prepare(double *a, int _n) {
    n = _n;
    nblocks = (n + B - 1) / B;
    _a = a;
    btree = (double(*)[8]) aligned_alloc(64, 64 * nblocks);
    build();
}

int cmp(dreg x_vec, double* y_ptr) {
    dreg y_vec = _mm256_load_pd(y_ptr);
    dreg mask = _mm256_cmp_pd(x_vec, y_vec, _CMP_GT_OQ);
    return _mm256_movemask_pd(mask);
}

double lower_bound(double x) {
    int k = 0;
    double res = INF;
    dreg x_vec = _mm256_set1_pd(x);
    while (k < nblocks) {
        int mask = ~(
            cmp(x_vec, &btree[k][0]) +
            (cmp(x_vec, &btree[k][4]) << 4)
        );
        int i = __builtin_ffs(mask) - 1;
        if (i < B)
            res = btree[k][i];
        k = go(k, i);
    }
    return res;
}

