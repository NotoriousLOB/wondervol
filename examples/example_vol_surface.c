#include <dispersion_scanner.h>
#include <stdio.h>
#include <stdlib.h>

const char* expiry_label(double T) {
    static char buf[32];
    if (T <= 1.0/12 + 0.01) snprintf(buf, sizeof(buf), "1M");
    else if (T <= 3.0/12 + 0.01) snprintf(buf, sizeof(buf), "3M");
    else if (T <= 6.0/12 + 0.01) snprintf(buf, sizeof(buf), "6M");
    else snprintf(buf, sizeof(buf), "1Y");
    return buf;
}

int main(void) {
    printf("=== Volatility Surface Example ===\n\n");
    
    double S = 100.0;
    double r = 0.03;
    double q = 0.0;
    
    // Strikes from 80 to 120
    double strikes[] = {80, 85, 90, 95, 100, 105, 110, 115, 120};
    size_t n_strikes = sizeof(strikes) / sizeof(strikes[0]);
    
    // Expiries from 1 month to 1 year
    double expiries[] = {1.0/12, 3.0/12, 6.0/12, 1.0};
    size_t n_expiries = sizeof(expiries) / sizeof(expiries[0]);
    
    // Base volatility
    double base_vol = 0.25;
    
    printf("Generating volatility surface for spot = $%.2f\n\n", S);
    
    // Simulate market prices with volatility smile
    printf("%-10s", "Strike");
    for (size_t t = 0; t < n_expiries; t++) {
        printf("%-12s", expiry_label(expiries[t]));
    }
    printf("\n");
    
    printf("%-10s", "------");
    for (size_t t = 0; t < n_expiries; t++) {
        printf("%-12s", "----");
    }
    printf("\n");
    
    for (size_t k = 0; k < n_strikes; k++) {
        printf("$%-9.0f", strikes[k]);
        
        for (size_t t = 0; t < n_expiries; t++) {
            // Generate market price with smile effect
            // Far OTM/ITM options have higher IV (volatility smile)
            double moneyness = strikes[k] / S;
            double smile_adj = 0.05 * (moneyness - 1.0) * (moneyness - 1.0);
            double time_adj = 0.02 * sqrt(expiries[t]);  // Term structure
            
            double true_vol = base_vol + smile_adj + time_adj;
            
            // Calculate option price with this vol
            double price = bsm_price(S, strikes[k], expiries[t], r, q, true_vol, 1);
            
            // Recover implied vol
            double iv = bsm_implied_vol(S, strikes[k], expiries[t], r, q, price, 1);
            
            printf("%-12.2f%%", iv * 100);
        }
        printf("\n");
    }
    
    printf("\n--- IV Recovery Accuracy ---\n");
    printf("Verifying implied vol calculation is accurate...\n");
    
    int n_tests = 0;
    int n_passed = 0;
    
    for (double strike = 70; strike <= 130; strike += 10) {
        for (double T = 0.25; T <= 2.0; T += 0.25) {
            for (double vol = 0.1; vol <= 0.5; vol += 0.1) {
                double price = bsm_price(S, strike, T, r, q, vol, 1);
                double iv = bsm_implied_vol(S, strike, T, r, q, price, 1);
                
                double error = iv - vol;
                if (fabs(error) < 1e-8) {
                    n_passed++;
                }
                n_tests++;
            }
        }
    }
    
    printf("Passed: %d/%d tests (%.1f%%)\n", n_passed, n_tests, 
           100.0 * n_passed / n_tests);
    
    printf("\n--- Vectorized IV Calculation ---\n");
    printf("Processing 10,000 options with vectorized solver...\n");
    
    const size_t n_opts = 10000;
    double* S_arr = aligned_alloc(64, n_opts * sizeof(double));
    double* K_arr = aligned_alloc(64, n_opts * sizeof(double));
    double* tau_arr = aligned_alloc(64, n_opts * sizeof(double));
    double* r_arr = aligned_alloc(64, n_opts * sizeof(double));
    double* q_arr = aligned_alloc(64, n_opts * sizeof(double));
    double* price_arr = aligned_alloc(64, n_opts * sizeof(double));
    double* iv_arr = aligned_alloc(64, n_opts * sizeof(double));
    uint8_t* is_call = aligned_alloc(64, n_opts * sizeof(uint8_t));
    
    // Generate random options
    for (size_t i = 0; i < n_opts; i++) {
        S_arr[i] = 100.0;
        K_arr[i] = 80.0 + (i % 50);
        tau_arr[i] = 0.1 + (i % 20) / 10.0;
        r_arr[i] = 0.03;
        q_arr[i] = 0.0;
        is_call[i] = (i % 2);
        price_arr[i] = bsm_price(S_arr[i], K_arr[i], tau_arr[i], r_arr[i], q_arr[i], 0.25, is_call[i]);
    }
    
    bsm_implied_vol_vec_arm(S_arr, K_arr, tau_arr, r_arr, q_arr, price_arr, is_call, n_opts, iv_arr);
    
    printf("Calculated %zu implied volatilities\n", n_opts);
    printf("Sample output (first 5):\n");
    for (size_t i = 0; i < 5 && i < n_opts; i++) {
        printf("  Option %zu: IV = %.4f%%\n", i+1, iv_arr[i] * 100);
    }
    
    free(S_arr); free(K_arr); free(tau_arr); free(r_arr); free(q_arr);
    free(price_arr); free(iv_arr); free(is_call);
    
    return 0;
}