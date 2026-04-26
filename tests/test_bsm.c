#include <wondervol.h>
#include <stdio.h>
#include <math.h>

#define TEST_EPSILON 1e-10

static int tests_passed = 0;
static int tests_failed = 0;

void test_assert(const char* name, int condition) {
    if (condition) {
        printf("✓ %s\n", name);
        tests_passed++;
    } else {
        printf("✗ %s\n", name);
        tests_failed++;
    }
}

void test_assert_double(const char* name, double actual, double expected, double epsilon) {
    int condition = fabs(actual - expected) < epsilon;
    if (condition) {
        printf("✓ %s (got %.10f, expected %.10f)\n", name, actual, expected);
        tests_passed++;
    } else {
        printf("✗ %s (got %.10f, expected %.10f, diff %.10f)\n", 
               name, actual, expected, fabs(actual - expected));
        tests_failed++;
    }
}

int main(void) {
    printf("=== BSM Pricing Tests ===\n\n");
    
    // Test 1: ATM call option — Cody erfc CDF has ~1.5e-15 absolute error → ~1e-12 price error
    double price_atm_call = bsm_price(100.0, 100.0, 1.0, 0.05, 0.0, 0.2, 1);
    test_assert_double("ATM call (S=K=100, T=1, r=5%, vol=20%)",
                       price_atm_call, 10.4505835721841546, 1e-10);

    // Test 2: ATM put option
    double price_atm_put = bsm_price(100.0, 100.0, 1.0, 0.05, 0.0, 0.2, 0);
    test_assert_double("ATM put (S=K=100, T=1, r=5%, vol=20%)",
                       price_atm_put, 5.5735260221841535, 1e-10);

    // Test 3: Deep ITM call
    double price_itm_call = bsm_price(150.0, 100.0, 1.0, 0.05, 0.0, 0.2, 1);
    double intrinsic_itm = 150.0 - 100.0 * exp(-0.05 * 1.0);
    test_assert("Deep ITM call > intrinsic", price_itm_call > intrinsic_itm);

    // Test 4: Deep OTM call (S=50, K=100 — still has some value with vol=20%)
    double price_otm_call = bsm_price(50.0, 100.0, 1.0, 0.05, 0.0, 0.2, 1);
    test_assert("Deep OTM call near zero", price_otm_call < 0.01 && price_otm_call >= 0.0);
    
    // Test 5: Put-call parity
    double S = 100.0, K = 95.0, T = 0.5, r = 0.03, q = 0.01, sigma = 0.25;
    double call = bsm_price(S, K, T, r, q, sigma, 1);
    double put = bsm_price(S, K, T, r, q, sigma, 0);
    double lhs = call + K * exp(-r * T);
    double rhs = put + S * exp(-q * T);
    test_assert_double("Put-call parity", lhs, rhs, 1e-10);
    
    // Test 6: Zero time to expiry - call
    double price_expired_call = bsm_price(110.0, 100.0, 0.0, 0.05, 0.0, 0.2, 1);
    test_assert_double("Expired call = intrinsic", price_expired_call, 10.0, 1e-10);
    
    // Test 7: Zero time to expiry - put
    double price_expired_put = bsm_price(90.0, 100.0, 0.0, 0.05, 0.0, 0.2, 0);
    test_assert_double("Expired put = intrinsic", price_expired_put, 10.0, 1e-10);
    
    // Test 8: Zero volatility
    double price_zero_vol = bsm_price(100.0, 100.0, 1.0, 0.05, 0.0, 0.0, 1);
    double fwd_zero_vol = 100.0 * exp(-0.0 * 1.0) - 100.0 * exp(-0.05 * 1.0);
    test_assert_double("Zero vol call = discounted fwd", price_zero_vol, fwd_zero_vol, 1e-10);
    
    // Test 9: With dividend yield
    double price_div_call = bsm_price(100.0, 100.0, 1.0, 0.05, 0.03, 0.2, 1);
    test_assert("Dividend reduces call price", price_div_call < price_atm_call);
    
    // Test 10: With dividend - put increases
    double price_div_put = bsm_price(100.0, 100.0, 1.0, 0.05, 0.03, 0.2, 0);
    double price_no_div_put = bsm_price(100.0, 100.0, 1.0, 0.05, 0.0, 0.2, 0);
    test_assert("Dividend increases put price", price_div_put > price_no_div_put);

    // === Edge Cases (C1) ===
    printf("\n--- Edge Cases ---\n");

    // Test 11: Extreme S/K ratio (deep ITM call)
    double price_extreme_itm = bsm_price(1000.0, 1.0, 1.0, 0.05, 0.0, 0.2, 1);
    test_assert("Extreme ITM call (S=1000, K=1) > spot intrinsic",
                price_extreme_itm > fmax(1000.0 - 1.0, 0.0));

    // Test 12: Extreme S/K ratio (deep OTM call)
    double price_extreme_otm = bsm_price(1.0, 1000.0, 1.0, 0.05, 0.0, 0.2, 1);
    test_assert("Extreme OTM call (S=1, K=1000) near zero",
                price_extreme_otm >= 0.0 && price_extreme_otm < 1e-10);

    // Test 13: Negative interest rate
    double price_neg_r = bsm_price(100.0, 100.0, 1.0, -0.02, 0.0, 0.2, 1);
    double price_pos_r = bsm_price(100.0, 100.0, 1.0, 0.02, 0.0, 0.2, 1);
    test_assert("Negative rate lowers call price", price_neg_r < price_pos_r);

    // Test 14: Put-call parity with negative rate
    {
        double Sn = 100.0, Kn = 100.0, Tn = 1.0, rn = -0.02, qn = 0.0, sign = 0.3;
        double cn = bsm_price(Sn, Kn, Tn, rn, qn, sign, 1);
        double pn = bsm_price(Sn, Kn, Tn, rn, qn, sign, 0);
        double lhsn = cn + Kn * exp(-rn * Tn);
        double rhsn = pn + Sn * exp(-qn * Tn);
        test_assert_double("Put-call parity (r=-2%)", lhsn, rhsn, 1e-10);
    }

    // Test 15: Very small T (near-expiry)
    double price_tiny_t = bsm_price(101.0, 100.0, 1e-6, 0.05, 0.0, 0.2, 1);
    test_assert_double("Near-expiry ITM call ≈ intrinsic", price_tiny_t, 1.0, 1e-3);

    // Test 16: Very high vol
    double price_high_vol = bsm_price(100.0, 100.0, 1.0, 0.05, 0.0, 5.0, 1);
    test_assert("High vol (σ=5.0) call is positive", price_high_vol > 0.0);
    test_assert("High vol call < spot", price_high_vol <= 100.0);

    // === CDF Accuracy (Cody Algorithm 715) ===
    printf("\n--- CDF Accuracy Tests (Cody Algorithm 715) ---\n");
    {
        /* High-precision reference values from Wolfram Alpha / scipy.special.ndtr */
        struct { double x; double expected; } cases[] = {
            { 0.0,   0.5                   },
            { 1.0,   0.8413447460685429    },
            {-1.0,   0.15865525393145705   },
            { 3.0,   0.9986501019683699    },
            {-3.0,   0.0013498980316301035 },
            { 6.0,   0.999999999013412     },
            {-6.0,   9.86587645037702e-10  },
        };
        int n_cases = (int)(sizeof(cases) / sizeof(cases[0]));
        for (int i = 0; i < n_cases; i++) {
            double got = norm_cdf(cases[i].x);
            char label[64];
            snprintf(label, sizeof(label), "norm_cdf(%.1f)", cases[i].x);
            test_assert_double(label, got, cases[i].expected, 1e-13);
        }
    }

    // === NEON Price Cross-validation ===
    printf("\n--- NEON vs Scalar Price Cross-validation ---\n");
#if defined(WV_ARCH_NEON)
    {
        // ATM
        double S_v[2] = {100.0, 100.0};
        double K_v[2] = {100.0, 100.0};
        double tau_v[2] = {1.0, 1.0};
        double r_v[2] = {0.05, 0.05};
        double q_v[2] = {0.0, 0.0};
        double sig_v[2] = {0.2, 0.2};

        float64x2_t neon_price = bsm_price_neon(
            vld1q_f64(S_v), vld1q_f64(K_v), vld1q_f64(tau_v),
            vld1q_f64(r_v), vld1q_f64(q_v), vld1q_f64(sig_v),
            vdupq_n_u64(~0ULL));  /* both calls */
        double np0 = vgetq_lane_f64(neon_price, 0);
        double sp0 = bsm_price(100.0, 100.0, 1.0, 0.05, 0.0, 0.2, 1);
        test_assert_double("NEON vs scalar: ATM call", np0, sp0, 1e-12);

        // Put
        float64x2_t neon_put = bsm_price_neon(
            vld1q_f64(S_v), vld1q_f64(K_v), vld1q_f64(tau_v),
            vld1q_f64(r_v), vld1q_f64(q_v), vld1q_f64(sig_v),
            vdupq_n_u64(0));  /* both puts */
        double nput0 = vgetq_lane_f64(neon_put, 0);
        double sput0 = bsm_price(100.0, 100.0, 1.0, 0.05, 0.0, 0.2, 0);
        test_assert_double("NEON vs scalar: ATM put", nput0, sput0, 1e-12);
    }
    {
        // ITM call + OTM put (mixed lanes)
        double S_v[2] = {150.0, 50.0};
        double K_v[2] = {100.0, 100.0};
        double tau_v[2] = {1.0, 0.5};
        double r_v[2] = {0.05, 0.03};
        double q_v[2] = {0.02, 0.0};
        double sig_v[2] = {0.3, 0.25};

        uint64x2_t mask = vsetq_lane_u64(~0ULL, vdupq_n_u64(0), 0); /* lane0=call, lane1=put */
        float64x2_t neon_pr = bsm_price_neon(
            vld1q_f64(S_v), vld1q_f64(K_v), vld1q_f64(tau_v),
            vld1q_f64(r_v), vld1q_f64(q_v), vld1q_f64(sig_v), mask);

        double np0 = vgetq_lane_f64(neon_pr, 0);
        double np1 = vgetq_lane_f64(neon_pr, 1);
        double sp0 = bsm_price(150.0, 100.0, 1.0, 0.05, 0.02, 0.3, 1);
        double sp1 = bsm_price(50.0, 100.0, 0.5, 0.03, 0.0, 0.25, 0);
        test_assert_double("NEON vs scalar: ITM call (mixed)", np0, sp0, 1e-12);
        test_assert_double("NEON vs scalar: OTM put (mixed)", np1, sp1, 1e-12);
    }
    {
        // Expired options
        double S_v[2] = {110.0, 90.0};
        double K_v[2] = {100.0, 100.0};
        double tau_v[2] = {0.0, 0.0};
        double r_v[2] = {0.05, 0.05};
        double q_v[2] = {0.0, 0.0};
        double sig_v[2] = {0.2, 0.2};

        uint64x2_t mask = vsetq_lane_u64(~0ULL, vdupq_n_u64(~0ULL), 1);
        float64x2_t neon_exp = bsm_price_neon(
            vld1q_f64(S_v), vld1q_f64(K_v), vld1q_f64(tau_v),
            vld1q_f64(r_v), vld1q_f64(q_v), vld1q_f64(sig_v), mask);

        test_assert_double("NEON expired: ITM call = 10", vgetq_lane_f64(neon_exp, 0), 10.0, 1e-12);
        test_assert_double("NEON expired: OTM call = 0", vgetq_lane_f64(neon_exp, 1), 0.0, 1e-12);
    }
#else
    printf("  (NEON not available on this platform — skipping NEON cross-validation)\n");
#endif

    printf("\n=== Summary ===\n");
    printf("Passed: %d\n", tests_passed);
    printf("Failed: %d\n", tests_failed);
    
    return tests_failed > 0 ? 1 : 0;
}
