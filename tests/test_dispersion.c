#include <wondervol.h>
#include <stdio.h>
#include <math.h>

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
    printf("=== Dispersion Scanner Tests ===\n\n");
    
    // Test 1: Single component correlation = 1.0
    ComponentVol single_comp = {1.0, 0.2, {0}};
    double rho1 = implied_correlation_arm(0.2, 1, &single_comp);
    test_assert_double("Single component: rho = 1.0", rho1, 1.0, 1e-10);
    
    // Test 2: Two identical components with same vol
    ComponentVol two_comps[2] = {
        {0.5, 0.2, {0}},
        {0.5, 0.2, {0}}
    };
    double etf_vol_same = 0.2;
    double rho2 = implied_correlation_arm(etf_vol_same, 2, two_comps);
    test_assert_double("Two identical components", rho2, 1.0, 1e-10);
    
    // Test 3: Zero correlation case
    // If ETF vol equals weighted RMS vol, correlation = 0
    ComponentVol uncorr_comps[2] = {
        {0.5, 0.3, {0}},
        {0.5, 0.1, {0}}
    };
    double sum_w2_sig2 = 0.5*0.5*0.3*0.3 + 0.5*0.5*0.1*0.1;
    double etf_vol_uncorr = sqrt(sum_w2_sig2);
    double rho3 = implied_correlation_arm(etf_vol_uncorr, 2, uncorr_comps);
    test_assert_double("Zero correlation case", rho3, 0.0, 1e-10);
    
    // Test 4: Correlation clamped to [-1, 1]
    double rho_high = implied_correlation_arm(0.5, 2, two_comps);
    test_assert("High correlation clamped to <= 1.0", rho_high <= 1.0);
    test_assert("High correlation clamped to >= -1.0", rho_high >= -1.0);
    
    // Test 5: Empty components returns NAN
    double rho_empty = implied_correlation_arm(0.2, 0, two_comps);
    test_assert("Empty components returns NAN", isnan(rho_empty));
    
    // Test 6: Dispersion opportunity scanning
    printf("\n--- Dispersion Opportunity Tests ---\n");
    
    size_t n_etf = 5;
    double etf_ivs[5] = {0.15, 0.18, 0.20, 0.22, 0.25};
    double etf_vegas[5] = {0.1, 0.12, 0.15, 0.18, 0.2};
    
    ComponentVol basket[3] = {
        {0.4, 0.25, {0}},
        {0.35, 0.22, {0}},
        {0.25, 0.20, {0}}
    };
    
    DispersionOpportunity opps[5];
    double target_rho = 0.5;
    double thresh = 0.1;
    double basket_vega = 1.0;
    
    scan_dispersion_opps_arm(etf_ivs, etf_vegas, n_etf, basket, 3, 
                             target_rho, thresh, basket_vega, opps);
    
    test_assert("Opportunities computed", opps[0].score > 0);
    
    // Find the best opportunity
    int best_idx = 0;
    double best_score = 0;
    for (size_t i = 0; i < n_etf; i++) {
        if (opps[i].score > best_score) {
            best_score = opps[i].score;
            best_idx = (int)i;
        }
    }
    printf("Best opportunity: ETF %d with score %.4f\n", best_idx, best_score);
    
    // Test 7: Single component basket
    DispersionOpportunity single_opp[2];
    scan_dispersion_opps_arm(etf_ivs, etf_vegas, 2, basket, 1, 
                             target_rho, thresh, basket_vega, single_opp);
    test_assert_double("Single component basket: rho = 1.0", 
                       single_opp[0].implied_rho, 1.0, 1e-10);
    
    // Test 8: Threshold filtering
    int above_thresh = 0;
    for (size_t i = 0; i < n_etf; i++) {
        if (opps[i].score > thresh) {
            above_thresh++;
        }
    }
    printf("Opportunities above threshold: %d\n", above_thresh);
    
    // Test 9: Vega neutrality calculation
    test_assert("Suggested notional computed", opps[0].suggested_etf_notional != 0.0);
    test_assert("Basket vega stored", opps[0].basket_vega == basket_vega);
    
    // Test 10: Multiple component basket correlation consistency
    double manual_rho = implied_correlation_arm(etf_ivs[2], 3, basket);
    test_assert_double("Manual correlation matches", opps[2].implied_rho, manual_rho, 1e-10);
    
    printf("\n=== Summary ===\n");
    printf("Passed: %d\n", tests_passed);
    printf("Failed: %d\n", tests_failed);
    
    return tests_failed > 0 ? 1 : 0;
}
