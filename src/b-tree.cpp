
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

void prepare(int *a, int n);
int lower_bound(int x);

int main(int argc, char* argv[]) {
    int n = (argc > 1 ? atoi(argv[1]) : N);
    int m = (argc > 2 ? atoi(argv[2]) : 1<<20);

    int *a = new int[n];
    int *q = new int[m];

    for (int i = 0; i < n; i++)
        a[i] = rand();
    for (int i = 0; i < m; i++)
        q[i] = rand();

    a[0] = RAND_MAX;
    std::sort(a, a + n);

    prepare(a, n);

    int checksum = 0;
    clock_t start = clock();

    #ifdef LATENCY
    int last = 0;

    for (int i = 0; i < m; i++) {
        last = lower_bound(q[i] ^ last);
        checksum ^= last;
    }
    #else
    for (int i = 0; i < m; i++)
        checksum ^= lower_bound(q[i]);
    #endif

    float seconds = float(clock() - start) / CLOCKS_PER_SEC;

    printf("%.2f ns per query\n", 1e9 * seconds / m);
    printf("%d\n", checksum);

    return 0;
}

#endif

typedef __m256i reg;

const int B = 16;
const int INF = std::numeric_limits<int>::max();

int n;
int nblocks;
int *_a;
int (*btree)[B];

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

void prepare(int *a, int _n) {
    n = _n;
    nblocks = (n + B - 1) / B;
    _a = a;
    btree = (int(*)[16]) aligned_alloc(64, 64 * nblocks);
    build();
}

int cmp(reg x_vec, int* y_ptr) {
    reg y_vec = _mm256_load_si256((reg*) y_ptr);
    reg mask = _mm256_cmpgt_epi32(x_vec, y_vec);
    return _mm256_movemask_ps((__m256) mask);
}

int lower_bound(int x) {
    int k = 0, res = INF;
    reg x_vec = _mm256_set1_epi32(x);
    while (k < nblocks) {
        int mask = ~(
            cmp(x_vec, &btree[k][0]) +
            (cmp(x_vec, &btree[k][8]) << 8)
        );
        int i = __builtin_ffs(mask) - 1;
        if (i < B)
            res = btree[k][i];
        k = go(k, i);
    }
    return res;
}

