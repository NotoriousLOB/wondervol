#include <dispersion_scanner.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#ifdef __APPLE__
#include <mach/mach_time.h>
static double get_time_ns(void) {
    static mach_timebase_info_data_t timebase;
    if (timebase.denom == 0) mach_timebase_info(&timebase);
    return (double)mach_absolute_time() * timebase.numer / timebase.denom;
}
#else
static double get_time_ns(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec * 1e9 + ts.tv_nsec;
}
#endif

static double random_double(double min, double max) {
    return min + (max - min) * ((double)rand() / RAND_MAX);
}

// Scalar IV solver — single-threaded baseline
static void bsm_implied_vol_scalar_st(const double* S, const double* K, const double* tau,
                                      const double* r, const double* q, const double* mkt_price,
                                      const uint8_t* is_call, size_t n, double* iv_out) {
    for (size_t i = 0; i < n; i++) {
        iv_out[i] = bsm_implied_vol(S[i], K[i], tau[i], r[i], q[i], mkt_price[i], is_call[i]);
    }
}

// Scalar IV solver — with OpenMP
static void bsm_implied_vol_scalar_omp(const double* S, const double* K, const double* tau,
                                       const double* r, const double* q, const double* mkt_price,
                                       const uint8_t* is_call, size_t n, double* iv_out) {
    #ifdef _OPENMP
    #pragma omp parallel for schedule(static, 256)
    #endif
    for (size_t i = 0; i < n; i++) {
        iv_out[i] = bsm_implied_vol(S[i], K[i], tau[i], r[i], q[i], mkt_price[i], is_call[i]);
    }
}

int main(int argc, char** argv) {
    size_t n = (argc > 1) ? atoi(argv[1]) : 100000;
    int iterations = (argc > 2) ? atoi(argv[2]) : 10;
    
    printf("=== NEON vs Scalar Benchmark ===\n");
    printf("Options count: %zu, Iterations: %d\n\n", n, iterations);
    
    // Allocate aligned arrays
    double* S = aligned_alloc(64, n * sizeof(double));
    double* K = aligned_alloc(64, n * sizeof(double));
    double* tau = aligned_alloc(64, n * sizeof(double));
    double* r = aligned_alloc(64, n * sizeof(double));
    double* q = aligned_alloc(64, n * sizeof(double));
    double* prices = aligned_alloc(64, n * sizeof(double));
    double* iv_neon = aligned_alloc(64, n * sizeof(double));
    double* iv_scalar = aligned_alloc(64, n * sizeof(double));
    uint8_t* is_call = aligned_alloc(64, n * sizeof(uint8_t));
    
    // Generate random options data
    srand(42);
    for (size_t i = 0; i < n; i++) {
        S[i] = random_double(50.0, 150.0);
        K[i] = random_double(80.0, 120.0);
        tau[i] = random_double(0.1, 2.0);
        r[i] = random_double(0.01, 0.08);
        q[i] = random_double(0.0, 0.05);
        is_call[i] = rand() % 2;
        
        double sigma = random_double(0.1, 0.5);
        prices[i] = bsm_price(S[i], K[i], tau[i], r[i], q[i], sigma, is_call[i]);
    }
    
    printf("Running %d iterations for each implementation...\n\n", iterations);
    
    // Warmup
    bsm_implied_vol_vec_arm(S, K, tau, r, q, prices, is_call, n, iv_neon);
    bsm_implied_vol_scalar_st(S, K, tau, r, q, prices, is_call, n, iv_scalar);

    // Benchmark 1: Scalar single-threaded
    double scalar_st_total = 0.0;
    for (int iter = 0; iter < iterations; iter++) {
        double start = get_time_ns();
        bsm_implied_vol_scalar_st(S, K, tau, r, q, prices, is_call, n, iv_scalar);
        double end = get_time_ns();
        scalar_st_total += (end - start) / 1e6;
    }
    double scalar_st_avg = scalar_st_total / iterations;

    // Benchmark 2: Scalar + OpenMP
    double scalar_omp_total = 0.0;
    for (int iter = 0; iter < iterations; iter++) {
        double start = get_time_ns();
        bsm_implied_vol_scalar_omp(S, K, tau, r, q, prices, is_call, n, iv_scalar);
        double end = get_time_ns();
        scalar_omp_total += (end - start) / 1e6;
    }
    double scalar_omp_avg = scalar_omp_total / iterations;

    // Benchmark 3: NEON + OpenMP
    double neon_total = 0.0;
    for (int iter = 0; iter < iterations; iter++) {
        double start = get_time_ns();
        bsm_implied_vol_vec_arm(S, K, tau, r, q, prices, is_call, n, iv_neon);
        double end = get_time_ns();
        neon_total += (end - start) / 1e6;
    }
    double neon_avg = neon_total / iterations;

    printf("=== Results ===\n");
    printf("Scalar (ST)   : %.3f ms (%.0f opts/sec)\n", scalar_st_avg, n / (scalar_st_avg / 1000.0));
    printf("Scalar (OMP)  : %.3f ms (%.0f opts/sec)\n", scalar_omp_avg, n / (scalar_omp_avg / 1000.0));
    printf("NEON + OMP    : %.3f ms (%.0f opts/sec)\n", neon_avg, n / (neon_avg / 1000.0));
    printf("\nSpeedup (NEON vs Scalar ST)  : %.2fx\n", scalar_st_avg / neon_avg);
    printf("Speedup (NEON vs Scalar OMP) : %.2fx\n", scalar_omp_avg / neon_avg);
    printf("Speedup (OMP vs ST)          : %.2fx\n", scalar_st_avg / scalar_omp_avg);

    // Verify NEON results via repricing
    printf("\n=== Correctness Check (NEON, via reprice) ===\n");
    double max_reprice_err = 0.0;
    for (size_t i = 0; i < n; i++) {
        double reprice = bsm_price(S[i], K[i], tau[i], r[i], q[i], iv_neon[i], is_call[i]);
        double err = fabs(reprice - prices[i]);
        if (err > max_reprice_err) max_reprice_err = err;
    }
    printf("Max reprice error: %.2e\n", max_reprice_err);
    printf("Accurate: %s\n", max_reprice_err < 1e-8 ? "YES" : "NO");
    
    // Cleanup
    free(S); free(K); free(tau); free(r); free(q);
    free(prices); free(iv_neon); free(iv_scalar); free(is_call);
    
    return 0;
}