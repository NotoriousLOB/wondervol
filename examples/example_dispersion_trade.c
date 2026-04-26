#include <wondervol.h>
#include <stdio.h>

int main(void) {
    printf("=== Dispersion Trading Example ===\n\n");
    
    // Define a basket of stocks (e.g., a sector ETF's components)
    ComponentVol basket[] = {
        {0.20, 0.28, {0}},  // Stock A: 20% weight, 28% vol
        {0.18, 0.25, {0}},  // Stock B: 18% weight, 25% vol
        {0.15, 0.22, {0}},  // Stock C: 15% weight, 22% vol
        {0.12, 0.30, {0}},  // Stock D: 12% weight, 30% vol
        {0.10, 0.24, {0}},  // Stock E: 10% weight, 24% vol
        {0.10, 0.26, {0}},  // Stock F: 10% weight, 26% vol
        {0.08, 0.32, {0}},  // Stock G: 8% weight, 32% vol
        {0.07, 0.35, {0}},  // Stock H: 7% weight, 35% vol
    };
    size_t n_comps = sizeof(basket) / sizeof(basket[0]);
    
    printf("Basket Components:\n");
    printf("%-10s %-10s %-10s\n", "Stock", "Weight", "IV");
    printf("%-10s %-10s %-10s\n", "-----", "------", "--");
    for (size_t i = 0; i < n_comps; i++) {
        printf("%-10zu %-10.2f %-10.2f%%\n", i+1, basket[i].weight, basket[i].iv * 100);
    }
    
    // Scenario: ETF is trading at different implied volatilities
    double etf_iv_scenarios[] = {0.15, 0.20, 0.22, 0.25, 0.30};
    size_t n_scenarios = sizeof(etf_iv_scenarios) / sizeof(etf_iv_scenarios[0]);
    
    printf("\n--- Implied Correlation Analysis ---\n");
    printf("%-15s %-15s %-15s\n", "ETF IV", "Implied Rho", "Opportunity");
    printf("%-15s %-15s %-15s\n", "------", "-----------", "-----------");
    
    double target_rho = 0.60;  // Expected correlation
    
    for (size_t i = 0; i < n_scenarios; i++) {
        double rho = implied_correlation_arm(etf_iv_scenarios[i], n_comps, basket);
        double deviation = rho - target_rho;
        const char* signal = (deviation > 0.1) ? "SHORT ETF" : 
                             (deviation < -0.1) ? "LONG ETF" : "NEUTRAL";
        
        printf("%-15.2f%% %-15.4f %-15s\n", etf_iv_scenarios[i] * 100, rho, signal);
    }
    
    // Full dispersion scan example
    printf("\n--- Full Dispersion Scan ---\n");
    
    // Multiple ETFs to scan
    double etf_ivs[] = {0.18, 0.22, 0.25, 0.28, 0.20};
    double etf_vegas[] = {0.12, 0.15, 0.18, 0.20, 0.14};
    size_t n_etfs = sizeof(etf_ivs) / sizeof(etf_ivs[0]);
    
    DispersionOpportunity opportunities[5];
    double threshold = 0.10;   // Min score to report
    double basket_vega = 2.0;  // Target basket vega exposure
    
    scan_dispersion_opps_arm(etf_ivs, etf_vegas, n_etfs, basket, n_comps,
                             target_rho, threshold, basket_vega, opportunities);
    
    printf("\nOpportunities (target rho = %.2f, threshold = %.2f):\n", 
           target_rho, threshold);
    printf("%-8s %-12s %-12s %-15s %-15s\n", 
           "ETF", "Implied Rho", "Score", "ETF Vega", "Suggested $");
    printf("%-8s %-12s %-12s %-15s %-15s\n", 
           "---", "-----------", "-----", "--------", "-----------");
    
    for (size_t i = 0; i < n_etfs; i++) {
        if (opportunities[i].score > 0) {
            printf("%-8zu %-12.4f %-12.4f %-15.4f $%-14.2f\n",
                   i+1,
                   opportunities[i].implied_rho,
                   opportunities[i].score,
                   opportunities[i].etf_vega,
                   opportunities[i].suggested_etf_notional);
        }
    }
    
    printf("\n--- Trade Strategy ---\n");
    printf("If implied correlation is HIGH (>%.2f):\n", target_rho);
    printf("  - Sell ETF options (high premium)\n");
    printf("  - Buy component options (cheaper)\n");
    printf("  - Expect correlation to revert lower\n");
    printf("\nIf implied correlation is LOW (<%.2f):\n", target_rho);
    printf("  - Buy ETF options (cheap)\n");
    printf("  - Sell component options (expensive)\n");
    printf("  - Expect correlation to revert higher\n");
    
    return 0;
}