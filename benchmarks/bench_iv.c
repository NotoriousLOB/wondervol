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

static int cmp_double(const void* a, const void* b) {
    double da = *(const double*)a, db = *(const double*)b;
    return (da > db) - (da < db);
}

static double random_double(double min, double max) {
    return min + (max - min) * ((double)rand() / RAND_MAX);
}

int main(int argc, char** argv) {
    size_t n = (argc > 1) ? atoi(argv[1]) : 100000;
    printf("=== Implied Volatility Benchmark ===\n");
    printf("Options count: %zu\n\n", n);
    
    // Allocate aligned arrays
    double* S = aligned_alloc(64, n * sizeof(double));
    double* K = aligned_alloc(64, n * sizeof(double));
    double* tau = aligned_alloc(64, n * sizeof(double));
    double* r = aligned_alloc(64, n * sizeof(double));
    double* q = aligned_alloc(64, n * sizeof(double));
    double* prices = aligned_alloc(64, n * sizeof(double));
    double* iv_out = aligned_alloc(64, n * sizeof(double));
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
    
    // Warmup
    printf("Warming up...\n");
    for (int i = 0; i < 3; i++) {
        bsm_implied_vol_vec_arm(S, K, tau, r, q, prices, is_call, n, iv_out);
    }
    
    // Benchmark
    printf("Running benchmark...\n");
    double start = get_time_ns();
    bsm_implied_vol_vec_arm(S, K, tau, r, q, prices, is_call, n, iv_out);
    double end = get_time_ns();
    
    double elapsed_ms = (end - start) / 1e6;
    double opts_per_sec = n / (elapsed_ms / 1000.0);
    
    printf("\n=== Results ===\n");
    printf("Time elapsed: %.3f ms\n", elapsed_ms);
    printf("Throughput: %.0f options/sec\n", opts_per_sec);
    printf("Time per option: %.3f μs\n", elapsed_ms * 1000.0 / n);
    
    // Verify accuracy across ALL samples
    printf("\n=== Accuracy Check (all %zu samples) ===\n", n);
    double max_error = 0.0;
    double total_error = 0.0;
    double p50_errors[101]; /* percentile buckets */
    for (int i = 0; i < 101; i++) p50_errors[i] = 0.0;

    double* all_errors = (double*)malloc(n * sizeof(double));

    for (size_t i = 0; i < n; i++) {
        double recalc_price = bsm_price(S[i], K[i], tau[i], r[i], q[i], iv_out[i], is_call[i]);
        double error = fabs(recalc_price - prices[i]);
        all_errors[i] = error;
        max_error = fmax(max_error, error);
        total_error += error;
    }

    qsort(all_errors, n, sizeof(double), cmp_double);

    printf("Max   price error: %.2e\n", max_error);
    printf("Mean  price error: %.2e\n", total_error / n);
    printf("p50   price error: %.2e\n", all_errors[n / 2]);
    printf("p90   price error: %.2e\n", all_errors[(size_t)(n * 0.9)]);
    printf("p99   price error: %.2e\n", all_errors[(size_t)(n * 0.99)]);
    printf("p99.9 price error: %.2e\n", all_errors[(size_t)(n * 0.999)]);

    free(all_errors);
    
    // Cleanup
    free(S); free(K); free(tau); free(r); free(q);
    free(prices); free(iv_out); free(is_call);
    
    return 0;
}