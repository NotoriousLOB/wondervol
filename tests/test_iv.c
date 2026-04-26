#include <wondervol.h>
#include <stdio.h>
#include <math.h>

#define TEST_EPSILON 1e-6

static int tests_passed = 0;
static int tests_failed = 0;

void test_assert_double(const char* name, double actual, double expected, double epsilon) {
    int condition = fabs(actual - expected) < epsilon || 
                    (isnan(expected) && isnan(actual));
    if (condition) {
        printf("✓ %s (got %.10f, expected %.10f)\n", name, actual, expected);
        tests_passed++;
    } else {
        printf("✗ %s (got %.10f, expected %.10f, diff %.10f)\n", 
               name, actual, expected, fabs(actual - expected));
        tests_failed++;
    }
}

void test_assert(const char* name, int condition) {
    if (condition) {
        printf("✓ %s\n", name);
        tests_passed++;
    } else {
        printf("✗ %s\n", name);
        tests_failed++;
    }
}

int main(void) {
    printf("=== Implied Volatility Tests ===\n\n");
    
    double S = 100.0, K = 100.0, T = 1.0, r = 0.05, q = 0.0;
    double true_vol = 0.25;
    
    // Test 1: Recover IV from BSM price
    double price = bsm_price(S, K, T, r, q, true_vol, 1);
    double iv = bsm_implied_vol(S, K, T, r, q, price, 1);
    test_assert_double("Recover IV from ATM call", iv, true_vol, 1e-10);

    // Test 2: Recover IV from put
    double put_price = bsm_price(S, K, T, r, q, true_vol, 0);
    double iv_put = bsm_implied_vol(S, K, T, r, q, put_price, 0);
    test_assert_double("Recover IV from ATM put", iv_put, true_vol, 1e-10);

    // Test 3: ITM option IV recovery
    double itm_price = bsm_price(120.0, 100.0, T, r, q, true_vol, 1);
    double iv_itm = bsm_implied_vol(120.0, 100.0, T, r, q, itm_price, 1);
    test_assert_double("Recover IV from ITM call", iv_itm, true_vol, 1e-10);

    // Test 4: OTM option IV recovery
    double otm_price = bsm_price(80.0, 100.0, T, r, q, true_vol, 1);
    double iv_otm = bsm_implied_vol(80.0, 100.0, T, r, q, otm_price, 1);
    test_assert_double("Recover IV from OTM call", iv_otm, true_vol, 1e-10);

    // Test 5: Short expiry
    double short_price = bsm_price(S, K, 0.1, r, q, true_vol, 1);
    double iv_short = bsm_implied_vol(S, K, 0.1, r, q, short_price, 1);
    test_assert_double("Recover IV from short expiry", iv_short, true_vol, 1e-9);
    
    // Test 6: Zero time to expiry
    double iv_zero = bsm_implied_vol(S, K, 0.0, r, q, 5.0, 1);
    test_assert_double("Zero expiry returns 0", iv_zero, 0.0, 1e-10);
    
    // Test 7: Vectorized IV solver with fixed vol
    printf("\n--- Vectorized IV Tests ---\n");
    size_t n = 1000;
    double S_arr[1000], K_arr[1000], tau_arr[1000], r_arr[1000], q_arr[1000];
    double price_arr[1000], iv_out[1000];
    uint8_t is_call[1000];

    for (size_t i = 0; i < n; i++) {
        S_arr[i] = 100.0 + (i % 20) * 2.5;
        K_arr[i] = 100.0;
        tau_arr[i] = 0.5 + (i % 10) * 0.1;
        r_arr[i] = 0.03 + (i % 5) * 0.01;
        q_arr[i] = 0.0;
        is_call[i] = (i % 2);
        price_arr[i] = bsm_price(S_arr[i], K_arr[i], tau_arr[i], r_arr[i], q_arr[i], 0.2, is_call[i]);
    }

    bsm_implied_vol_vec(S_arr, K_arr, tau_arr, r_arr, q_arr, price_arr, is_call, n, iv_out);

    double max_err = 0.0, sum_err = 0.0;
    for (size_t i = 0; i < n; i++) {
        double err = fabs(iv_out[i] - 0.2);
        if (err > max_err) max_err = err;
        sum_err += err;
    }
    printf("  n=1000 fixed vol: max_err=%.2e mean_err=%.2e\n", max_err, sum_err / n);
    test_assert("Vectorized IV recovery (1000 options, max_err < 1e-8)", max_err < 1e-8);

    // Test 8: Odd count (scalar tail)
    bsm_implied_vol_vec(S_arr, K_arr, tau_arr, r_arr, q_arr, price_arr, is_call, 999, iv_out);
    max_err = 0.0;
    for (size_t i = 0; i < 999; i++) {
        double err = fabs(iv_out[i] - 0.2);
        if (err > max_err) max_err = err;
    }
    test_assert("Vectorized IV with odd count (999 options)", max_err < 1e-8);

    // Test 9: n=1 (scalar tail only, no NEON iteration)
    {
        double S1[1] = {100.0}, K1[1] = {100.0}, tau1[1] = {1.0};
        double r1[1] = {0.05}, q1[1] = {0.0}, p1[1], iv1[1];
        uint8_t c1[1] = {1};
        p1[0] = bsm_price(100.0, 100.0, 1.0, 0.05, 0.0, 0.25, 1);
        bsm_implied_vol_vec(S1, K1, tau1, r1, q1, p1, c1, 1, iv1);
        test_assert_double("n=1 scalar tail", iv1[0], 0.25, 1e-10);
    }

    // Test 10: n=2 (one NEON iteration, no tail)
    {
        double S2[2] = {100.0, 110.0}, K2[2] = {100.0, 100.0}, tau2[2] = {1.0, 0.5};
        double r2[2] = {0.05, 0.03}, q2[2] = {0.0, 0.0}, p2[2], iv2[2];
        uint8_t c2[2] = {1, 0};
        p2[0] = bsm_price(100.0, 100.0, 1.0, 0.05, 0.0, 0.3, 1);
        p2[1] = bsm_price(110.0, 100.0, 0.5, 0.03, 0.0, 0.3, 0);
        bsm_implied_vol_vec(S2, K2, tau2, r2, q2, p2, c2, 2, iv2);
        test_assert_double("n=2 lane 0 (ATM call)", iv2[0], 0.3, 1e-10);
        test_assert_double("n=2 lane 1 (ITM put)", iv2[1], 0.3, 1e-10);
    }

    // Test 11: Varying true vols (0.05 to 2.0)
    printf("\n--- Varying Volatility Tests ---\n");
    {
        double test_vols[] = {0.05, 0.1, 0.15, 0.2, 0.3, 0.5, 0.8, 1.0, 1.5, 2.0};
        int n_vols = (int)(sizeof(test_vols) / sizeof(test_vols[0]));
        double Sv[2], Kv[2], tv[2], rv[2], qv[2], pv[2], ivv[2];
        uint8_t cv[2];

        for (int v = 0; v < n_vols; v++) {
            Sv[0] = 100.0; Sv[1] = 100.0;
            Kv[0] = 100.0; Kv[1] = 100.0;
            tv[0] = 1.0;   tv[1] = 1.0;
            rv[0] = 0.05;  rv[1] = 0.05;
            qv[0] = 0.0;   qv[1] = 0.0;
            cv[0] = 1;      cv[1] = 0;
            pv[0] = bsm_price(100.0, 100.0, 1.0, 0.05, 0.0, test_vols[v], 1);
            pv[1] = bsm_price(100.0, 100.0, 1.0, 0.05, 0.0, test_vols[v], 0);
            bsm_implied_vol_vec(Sv, Kv, tv, rv, qv, pv, cv, 2, ivv);

            char label[64];
            snprintf(label, sizeof(label), "Vol %.2f call recovery", test_vols[v]);
            test_assert_double(label, ivv[0], test_vols[v], 1e-8);
            snprintf(label, sizeof(label), "Vol %.2f put recovery", test_vols[v]);
            test_assert_double(label, ivv[1], test_vols[v], 1e-8);
        }
    }

    // Test 12: Deep OTM stress test
    printf("\n--- Deep OTM Stress Test ---\n");
    {
        double Sv[2] = {135.0, 145.0};
        double Kv[2] = {100.0, 100.0};
        double tv[2] = {0.9, 1.3};
        double rv[2] = {0.07, 0.06};
        double qv[2] = {0.0, 0.0};
        uint8_t cv[2] = {0, 0}; /* both deep OTM puts */
        double pv[2], ivv[2];
        pv[0] = bsm_price(135.0, 100.0, 0.9, 0.07, 0.0, 0.2, 0);
        pv[1] = bsm_price(145.0, 100.0, 1.3, 0.06, 0.0, 0.2, 0);
        bsm_implied_vol_vec(Sv, Kv, tv, rv, qv, pv, cv, 2, ivv);
        test_assert_double("Deep OTM put (S=135, K=100)", ivv[0], 0.2, 1e-8);
        test_assert_double("Deep OTM put (S=145, K=100)", ivv[1], 0.2, 1e-8);
    }

    // Test 13: Convergence failure — impossible prices (C3)
    printf("\n--- Convergence Failure Tests ---\n");
    {
        /* Call priced below intrinsic — solver should return clamped result, not crash */
        double iv_bad = bsm_implied_vol(110.0, 100.0, 1.0, 0.05, 0.0, 5.0, 1);
        test_assert("Below-intrinsic call: IV in [0.001, 5.0]",
                    iv_bad >= 0.001 && iv_bad <= 5.0);

        /* Negative price — solver should not crash */
        double iv_neg = bsm_implied_vol(100.0, 100.0, 1.0, 0.05, 0.0, -1.0, 1);
        test_assert("Negative price: IV in [0.001, 5.0]",
                    iv_neg >= 0.001 && iv_neg <= 5.0);

        /* Zero price */
        double iv_zero = bsm_implied_vol(100.0, 100.0, 1.0, 0.05, 0.0, 0.0, 1);
        test_assert("Zero price: IV in [0.001, 5.0]",
                    iv_zero >= 0.001 && iv_zero <= 5.0);
    }

    printf("\n=== Summary ===\n");
    printf("Passed: %d\n", tests_passed);
    printf("Failed: %d\n", tests_failed);

    return tests_failed > 0 ? 1 : 0;
}
