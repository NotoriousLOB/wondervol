#include <wondervol.h>
#include <stdio.h>
#include <math.h>

int main(void) {
    printf("=== Platform Detection Test ===\n\n");

#if defined(WV_ARCH_NEON)
    printf("SIMD backend : ARM NEON (2-wide float64x2_t)\n");
#elif defined(WV_ARCH_AVX512)
    printf("SIMD backend : AVX512 (8-wide __m512d)\n");
#elif defined(WV_ARCH_AVX2)
    printf("SIMD backend : AVX2 (4-wide __m256d)\n");
#else
    printf("SIMD backend : Scalar fallback\n");
#endif
    printf("\n");

    /* Cross-validate bsm_implied_vol_vec against scalar on 8 representative options */
    double S[8]   = {100, 110,  90, 105,  95, 100, 100, 100};
    double K[8]   = {100, 100, 100, 100, 100,  95, 105, 100};
    double tau[8] = {1.0, 1.0, 1.0, 0.5, 0.5, 1.0, 1.0, 0.25};
    double r[8]   = {0.05, 0.05, 0.05, 0.03, 0.03, 0.05, 0.05, 0.05};
    double q[8]   = {0, 0, 0, 0, 0, 0, 0, 0};
    uint8_t is_call[8] = {1, 1, 0, 1, 0, 1, 0, 1};
    double prices[8], iv_vec[8], iv_scalar[8];

    for (int i = 0; i < 8; i++)
        prices[i] = bsm_price(S[i], K[i], tau[i], r[i], q[i], 0.25, is_call[i]);

    bsm_implied_vol_vec(S, K, tau, r, q, prices, is_call, 8, iv_vec);

    for (int i = 0; i < 8; i++)
        iv_scalar[i] = bsm_implied_vol(S[i], K[i], tau[i], r[i], q[i], prices[i], is_call[i]);

    double max_diff = 0.0;
    for (int i = 0; i < 8; i++) {
        double d = fabs(iv_vec[i] - iv_scalar[i]);
        if (d > max_diff) max_diff = d;
    }
    printf("Vec vs scalar max diff: %.2e\n", max_diff);

    int pass = (max_diff < 1e-10);
    printf("%s vec vs scalar agreement (tol 1e-10)\n", pass ? "PASS" : "FAIL");
    printf("\n=== %s ===\n", pass ? "PASSED" : "FAILED");
    return pass ? 0 : 1;
}
