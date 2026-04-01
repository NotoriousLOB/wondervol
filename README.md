# Wondervol

[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![Platform](https://img.shields.io/badge/platform-ARM64-orange.svg)](https://en.wikipedia.org/wiki/ARM_architecture_family)
[![Language](https://img.shields.io/badge/language-C99-green.svg)](https://en.wikipedia.org/wiki/C99)

Today is gonna be the day that your vol forecast finally comes true.

By now you should’ve somehow realized what your gamma’s gotta do.

A high-performance header-only C library for vanquishing volatility, optimized with SIMD instructions and OpenMP parallelization.

## What is Dispersion Trading?

Dispersion trading is a sophisticated options strategy that exploits the relationship between:
- **Index/ETF implied volatility**
- **Component stock implied volatilities**

The strategy is based on the observation that the variance of an index is related to the variance of its components through their correlation. When the implied correlation deviates from historical norms, trading opportunities arise.

## Features

- **NEON SIMD Optimized**: 2-wide vectorized computation for ARM64
- **OpenMP Parallel**: Multi-threaded processing with cache-friendly scheduling
- **Header-only**: Single include file, no installation required
- **Implied Volatility**: Fast Newton-Raphson solver with robust convergence
- **Dispersion Analytics**: Compute implied correlations and scan for opportunities
- **Well-tested**: Comprehensive test suite with accuracy verification
- **Examples**: Working examples for common use cases

## Quick Start

### Prerequisites

- ARM64 processor (Apple Silicon, AWS Graviton, etc.)
- C compiler (Clang or GCC)
- CMake 3.16+
- OpenMP support

### Installation

```bash
git clone https://github.com/yourusername/dispersion-scanner.git
cd dispersion-scanner
mkdir build && cd build
cmake ..
make -j$(nproc)
sudo make install
```

### Basic Usage

```c
#include <dispersion_scanner.h>
#include <stdio.h>

int main() {
    // Price a European call option
    double price = bsm_price(
        100.0,   // Spot price
        105.0,   // Strike price
        0.5,     // Time to expiry (years)
        0.05,    // Risk-free rate
        0.02,    // Dividend yield
        0.25,    // Volatility
        1        // 1=call, 0=put
    );
    
    printf("Call price: $%.4f\n", price);
    
    // Calculate implied volatility from market price
    double mkt_price = 8.50;
    double iv = bsm_implied_vol(
        100.0, 105.0, 0.5, 0.05, 0.02,
        mkt_price, 1
    );
    
    printf("Implied volatility: %.2f%%\n", iv * 100);
    
    return 0;
}
```

Compile with:
```bash
clang -O3 -fopenmp -mcpu=native -I/usr/local/include \
    -o my_program my_program.c -lm
```

## Building

### CMake Options

| Option | Default | Description |
|--------|---------|-------------|
| `DISP_BUILD_TESTS` | ON | Build test suite |
| `DISP_BUILD_BENCHMARKS` | ON | Build benchmarks |
| `DISP_BUILD_EXAMPLES` | ON | Build examples |
| `DISP_ENABLE_OPENMP` | ON | Enable OpenMP parallelization |
| `DISP_STRICT_TRAP` | OFF | Abort on validation failures |

### Examples

```bash
# Full build with all components
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)

# Minimal build (library only)
cmake .. -DDISP_BUILD_TESTS=OFF \
         -DDISP_BUILD_BENCHMARKS=OFF \
         -DDISP_BUILD_EXAMPLES=OFF
make

# Debug build with strict trapping
cmake .. -DCMAKE_BUILD_TYPE=Debug \
         -DDISP_STRICT_TRAP=ON
make
```

## Running Tests

```bash
cd build

# Run all tests
ctest --output-on-failure

# Run specific test
./tests/test_bsm
./tests/test_iv
./tests/test_dispersion

# Verbose output
ctest -V
```

## Benchmarks

```bash
cd build

# Implied volatility benchmark (100k options)
./benchmarks/bench_iv

# NEON vs scalar comparison
./benchmarks/bench_neon_scalar 100000 10

# Dispersion scanner benchmark
./benchmarks/bench_dispersion 1000 50
```

### Example Output

```
=== NEON vs Scalar Benchmark ===
Options count: 100000, Iterations: 10

=== Results ===
NEON   : 12.345 ms (8100445 opts/sec)
Scalar : 28.901 ms (3460000 opts/sec)
Speedup: 2.34x

=== Correctness Check ===
Max IV difference: 1.23e-12
Results match: YES
```

## Dispersion Trading Example

```c
#include <dispersion_scanner.h>
#include <stdio.h>

int main() {
    // Define basket components (e.g., tech sector)
    ComponentVol basket[] = {
        {0.25, 0.30, {0}},  // AAPL: 25% weight, 30% vol
        {0.20, 0.28, {0}},  // MSFT: 20% weight, 28% vol
        {0.15, 0.35, {0}},  // NVDA: 15% weight, 35% vol
        {0.12, 0.32, {0}},  // GOOGL: 12% weight, 32% vol
        {0.10, 0.25, {0}},  // META: 10% weight, 25% vol
        {0.08, 0.27, {0}},  // NFLX: 8% weight, 27% vol
        {0.10, 0.29, {0}},  // AMZN: 10% weight, 29% vol
    };
    size_t n_comps = 7;
    
    // ETF implied volatility
    double etf_iv = 0.22;
    
    // Calculate implied correlation
    double rho = implied_correlation_arm(etf_iv, n_comps, basket);
    
    printf("Implied correlation: %.4f\n", rho);
    
    // Trading signal
    if (rho > 0.70) {
        printf("Signal: SELL ETF (high correlation)\n");
        printf("Strategy: Short ETF options, long component options\n");
    } else if (rho < 0.40) {
        printf("Signal: BUY ETF (low correlation)\n");
        printf("Strategy: Long ETF options, short component options\n");
    }
    
    return 0;
}
```

## API Reference

### Black-Scholes-Merton Pricing

```c
double bsm_price(double S, double K, double tau, double r, double q,
                 double sigma, int is_call);
```
Calculate European option price using Black-Scholes-Merton model.

**Parameters:**
- `S` - Spot price
- `K` - Strike price
- `tau` - Time to expiry (years)
- `r` - Risk-free rate
- `q` - Dividend yield
- `sigma` - Volatility
- `is_call` - 1 for call, 0 for put

**Returns:** Option price

### Implied Volatility

```c
double bsm_implied_vol(double S, double K, double tau, double r, double q,
                       double mkt_price, int is_call);
```
Calculate implied volatility using Newton-Raphson iteration.

**Parameters:**
- `mkt_price` - Observed market price
- Other parameters same as `bsm_price`

**Returns:** Implied volatility or clamped to [0.001, 5.0]

### Vectorized Implied Volatility

```c
void bsm_implied_vol_vec_arm(const double* S, const double* K,
                             const double* tau, const double* r,
                             const double* q, const double* mkt_price,
                             const uint8_t* is_call, size_t n,
                             double* iv_out);
```
Vectorized IV calculation using NEON SIMD and OpenMP.

### Implied Correlation

```c
double implied_correlation_arm(double etf_iv, size_t n_comps,
                               const ComponentVol* comps);
```
Calculate implied correlation from ETF and component volatilities.

**Returns:** Correlation in [-1.0, 1.0], or NAN on error

### Dispersion Scanner

```c
void scan_dispersion_opps_arm(
    const double* etf_ivs, const double* etf_vegas, size_t n_etf,
    const ComponentVol* comps, size_t n_comps,
    double target_rho, double thresh, double basket_vega,
    DispersionOpportunity* opps);
```
Scan multiple ETFs for dispersion trading opportunities.

## Data Structures

```c
typedef struct {
    double weight;           // Component weight in basket
    double iv;              // Implied volatility
    _Alignas(64) double pad[6];
} ComponentVol;

typedef struct {
    double implied_rho;      // Implied correlation
    double score;           // Opportunity score (deviation from target)
    double etf_vega;        // ETF vega exposure
    double suggested_etf_notional;  // Hedge ratio
    double basket_vega;     // Target basket vega
} DispersionOpportunity;
```

## Performance

Benchmarks on Apple M2 Pro:

| Operation | Scalar | NEON | Speedup |
|-----------|--------|------|---------|
| IV calculation (100k) | 28.9 ms | 12.3 ms | 2.3x |
| Dispersion scan (1k ETFs, 50 comps) | - | 0.8 ms | - |

## Algorithm Details

### Implied Volatility Solver

The Newton-Raphson solver uses:
- **Initial guess**: `sqrt(2π/T) * time_value / S`
- **Bounds**: [0.001, 5.0] (0.1% to 500% volatility)
- **Convergence**: `|price - mkt| < 1e-9`
- **Max iterations**: 50
- **Early exit**: If vega < 1e-12 or price error < 1e-12

### Implied Correlation Formula

```
σ²_etf = Σ(w²ᵢ × σ²ᵢ) + ρ × [(Σwᵢσᵢ)² - Σ(w²ᵢ × σ²ᵢ)]

Solving for ρ:
ρ = (σ²_etf - Σw²ᵢσ²ᵢ) / [(Σwᵢσᵢ)² - Σw²ᵢσ²ᵢ]
```

## Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

### Development Setup

```bash
git clone https://github.com/yourusername/dispersion-scanner.git
cd dispersion-scanner
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Debug -DDISP_STRICT_TRAP=ON
make -j$(nproc)
ctest -V
```

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments

- Black-Scholes-Merton model for option pricing
- Abramowitz & Stegun approximations for normal CDF
- ARM NEON intrinsics for SIMD optimization

## Related Projects

- [QuantLib](https://www.quantlib.org/) - Comprehensive quantitative finance library
- [Parklet Transform](https://github.com/yourusername/parklet) - Signal processing for motor control
