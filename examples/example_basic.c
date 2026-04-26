#include <wondervol.h>
#include <stdio.h>

int main(void) {
    printf("=== Dispersion Scanner - Basic Example ===\n\n");
    
    // Example 1: Price a single option
    double S = 100.0;      // Spot price
    double K = 105.0;      // Strike price
    double T = 0.5;        // Time to expiry (years)
    double r = 0.05;       // Risk-free rate
    double q = 0.02;       // Dividend yield
    double sigma = 0.25;   // Volatility
    
    double call_price = bsm_price(S, K, T, r, q, sigma, 1);
    double put_price = bsm_price(S, K, T, r, q, sigma, 0);
    
    printf("Option Pricing (Black-Scholes-Merton)\n");
    printf("  Spot:        $%.2f\n", S);
    printf("  Strike:      $%.2f\n", K);
    printf("  Time:        %.2f years\n", T);
    printf("  Rate:        %.1f%%\n", r * 100);
    printf("  Dividend:    %.1f%%\n", q * 100);
    printf("  Volatility:  %.1f%%\n", sigma * 100);
    printf("\n");
    printf("  Call price:  $%.4f\n", call_price);
    printf("  Put price:   $%.4f\n", put_price);
    
    // Example 2: Calculate implied volatility
    printf("\n--- Implied Volatility ---\n");
    double mkt_price = 8.50;  // Market observed price
    double iv = bsm_implied_vol(S, K, T, r, q, mkt_price, 1);
    printf("  Market price: $%.2f\n", mkt_price);
    printf("  Implied vol:  %.2f%%\n", iv * 100);
    
    // Example 3: Verify - price with implied vol should match market price
    double verify_price = bsm_price(S, K, T, r, q, iv, 1);
    printf("  Verify price: $%.4f (error: %.6f)\n", verify_price, 
           verify_price - mkt_price);
    
    return 0;
}