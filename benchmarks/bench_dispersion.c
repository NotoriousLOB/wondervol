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

int main(int argc, char** argv) {
    size_t n_etfs = (argc > 1) ? atoi(argv[1]) : 1000;
    size_t n_comps = (argc > 2) ? atoi(argv[2]) : 50;
    
    printf("=== Dispersion Scanner Benchmark ===\n");
    printf("ETFs: %zu, Components per basket: %zu\n\n", n_etfs, n_comps);
    
    // Allocate data
    double* etf_ivs = aligned_alloc(64, n_etfs * sizeof(double));
    double* etf_vegas = aligned_alloc(64, n_etfs * sizeof(double));
    DispersionOpportunity* opps = aligned_alloc(64, n_etfs * sizeof(DispersionOpportunity));
    
    ComponentVol* comps = aligned_alloc(64, n_comps * sizeof(ComponentVol));
    
    // Generate synthetic data
    srand(42);
    for (size_t i = 0; i < n_etfs; i++) {
        etf_ivs[i] = 0.15 + (rand() % 200) / 1000.0;  // 0.15 - 0.35
        etf_vegas[i] = 0.05 + (rand() % 100) / 1000.0;  // 0.05 - 0.15
    }
    
    double total_weight = 0.0;
    for (size_t i = 0; i < n_comps; i++) {
        comps[i].weight = 1.0 / n_comps;
        comps[i].iv = 0.15 + (rand() % 250) / 1000.0;
        total_weight += comps[i].weight;
    }
    
    // Normalize weights
    for (size_t i = 0; i < n_comps; i++) {
        comps[i].weight /= total_weight;
    }
    
    // Warmup
    printf("Warming up...\n");
    for (int i = 0; i < 5; i++) {
        scan_dispersion_opps_arm(etf_ivs, etf_vegas, n_etfs, comps, n_comps,
                                 0.5, 0.1, 1.0, opps);
    }
    
    // Benchmark
    printf("Running benchmark...\n");
    double start = get_time_ns();
    scan_dispersion_opps_arm(etf_ivs, etf_vegas, n_etfs, comps, n_comps,
                             0.5, 0.1, 1.0, opps);
    double end = get_time_ns();
    
    double elapsed_ms = (end - start) / 1e6;
    double scans_per_sec = n_etfs / (elapsed_ms / 1000.0);
    
    printf("\n=== Results ===\n");
    printf("Time elapsed: %.3f ms\n", elapsed_ms);
    printf("Throughput: %.0f ETFs/sec\n", scans_per_sec);
    printf("Time per ETF: %.3f μs\n", elapsed_ms * 1000.0 / n_etfs);
    
    // Find best opportunity
    int best_idx = 0;
    double best_score = 0.0;
    for (size_t i = 0; i < n_etfs; i++) {
        if (opps[i].score > best_score) {
            best_score = opps[i].score;
            best_idx = (int)i;
        }
    }
    
    printf("\nBest opportunity: ETF %d\n", best_idx);
    printf("  Implied correlation: %.4f\n", opps[best_idx].implied_rho);
    printf("  Score: %.4f\n", opps[best_idx].score);
    printf("  Suggested notional: %.2f\n", opps[best_idx].suggested_etf_notional);
    
    free(etf_ivs); free(etf_vegas); free(opps); free(comps);
    
    return 0;
}