# WonderVol

[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![Platform](https://img.shields.io/badge/platform-ARM64%20%7C%20x86--64%20%7C%20Scalar-orange.svg)](https://en.wikipedia.org/wiki/SIMD)
[![Language](https://img.shields.io/badge/language-C99-green.svg)](https://en.wikipedia.org/wiki/C99)

Today is gonna be the day that your vol forecast finally comes true.

By now you should've somehow realized what your gamma's gotta do.

Anyway, here's WonderVol: a high-performance header-only C library for vanquishing volatility and scanning dispersion trading opportunities, with platform-adaptive SIMD backends.

## What is Dispersion Trading?

Dispersion trading exploits the gap between index (ETF) implied volatility and component stock implied volatilities. The variance of a weighted basket equals:

```
σ²_etf = Σ(w²ᵢ σ²ᵢ) + ρ [(Σwᵢ σᵢ)² − Σ(w²ᵢ σ²ᵢ)]
```

When the implied correlation `ρ` extracted from market prices diverges from its historical or expected value, a dispersion trade captures the reversion:

- **High implied correlation**: sell ETF vol, buy single-stock vol
- **Low implied correlation**: buy ETF vol, sell single-stock vol

## Features

- **Cross-platform SIMD**: ARM NEON (2-wide), AVX2 (4-wide), AVX512 (8-wide), scalar fallback — detected at compile time, no runtime overhead
- **Header-only**: single `#include <wondervol.h>`, no compiled library required
- **Accurate normal CDF**: FDLIBM `precise_erfc` piecewise rational approximation, ~1 ulp error — replaces the A&S approximation which had ~7.5e-8 absolute error
- **Smart IV initial guess**: Corrado-Miller (1996) moneyness-adaptive formula, significantly reducing Newton-Raphson iterations for OTM options
- **Clean Newton-Raphson**: unified convergence threshold (1e-12), consistent vega-singularity handling across all backends, per-lane active masks in SIMD paths
- **OpenMP**: multi-threaded array IV solvers with static cache-friendly scheduling

## Mathematical Foundations

### Black-Scholes-Merton Pricing

European option price under the BSM model with continuous dividends:

```
d₁ = [ln(S/K) + (r − q + σ²/2)τ] / (σ√τ)
d₂ = d₁ − σ√τ

Call = S e^{-qτ} Φ(d₁)  − K e^{-rτ} Φ(d₂)
Put  = K e^{-rτ} Φ(−d₂) − S e^{-qτ} Φ(−d₁)
```

where `Φ` is the standard normal CDF, `S` spot, `K` strike, `τ` time to expiry, `r` risk-free rate, `q` continuous dividend yield, `σ` volatility.

### Normal CDF: FDLIBM `precise_erfc`

```
Φ(x) = erfc(−x/√2) / 2
```

`erfc` is computed via a piecewise rational approximation ported from the Sun/FreeBSD FDLIBM heritage, achieving ~1 ulp accuracy across all real inputs.
Four regions: `|x| < 0.84375`, `|x| < 1.25`, `|x| < 26.5`, and the asymptotic tail.
Each region uses a Horner-evaluated rational p/q form; the core of the `|x| < 26.5` paths uses `exp(−x² − 0.5625 + R(x)/S(x)) / x` to avoid catastrophic cancellation.

### Implied Volatility Solver

**Initial guess** — Corrado & Miller (1996), *J. Banking & Finance* 20:151–169:

```
C̃ = undiscounted call price (puts converted via parity: C̃ = P/e^{-rτ} + (F−K))
F = S e^{(r−q)τ}   (forward price)

σ₀ = √(2π/τ) / (S e^{-qτ}) × [C̃ − (F−K)/2 + √((C̃ − (F−K)/2)² − (F−K)²/π)]
```

Falls back to σ = 0.2 only when this formula yields a non-finite or negative result.

**Newton-Raphson iteration** on `f(σ) = BSM(σ) − C_mkt = 0`:

```
σₙ₊₁ = σₙ − f(σₙ) / f'(σₙ)
f'(σ) = vega = S e^{-qτ} Φ'(d₁) √τ
```

- Convergence: `|f(σ)| < 1e-12`
- Vega singularity stop: if `vega < 1e-12`, keep best estimate and exit
- Step damping: `|Δσ| ≤ σ/2` per iteration
- Bounds: `σ ∈ [0.001, 5.0]`
- Max iterations: 50 (scalar), 50 per lane (SIMD)

In SIMD paths, each lane carries its own active-mask; a lane is frozen when its vega drops below the threshold while other lanes continue iterating.

### Implied Correlation

Solving the basket-variance identity for `ρ`:

```
ρ = (σ²_etf − Σw²ᵢσ²ᵢ) / [(Σwᵢσᵢ)² − Σw²ᵢσ²ᵢ]
```

Result is clamped to `[−1, 1]`. Returns `NAN` for empty baskets; `1.0` for single-component baskets (perfect correlation by definition) and when the denominator is near zero.

## Limitations

The following are known limitations of the current implementation. We document them here so users don't repeat the mistakes of others who were foolish enough to deploy our code in production.

- **European options only.** BSM does not model early exercise; American options
  require different methods (binomial trees, finite differences).
- **Continuous dividend yield only.** Discrete dividends, dividend futures, and
  term-structure dividend models are not supported.
- **No smile or skew.** BSM assumes lognormal returns with constant volatility;
  real markets exhibit skew and excess kurtosis. Implied vols computed here are
  Black-Scholes flat-vol numbers, not stochastic vol parameters.
- **No term structure.** A single constant `σ` per option; no interpolation
  from a vol surface.
- **NEON/AVX exp and log are scalar lane-extract wrappers.** There is no native
  vectorized `exp`/`log` for `float64` in NEON or AVX2 without vendor libraries
  (ARM SVML, Intel SVML). The wrappers extract each lane to scalar, call
  `exp`/`log`, and repack. This is correct but slower than a polynomial
  approximation would be; the bottleneck is the Newton-Raphson CDF evaluation,
  not these calls.
- **Compile-time SIMD dispatch only.** Platform selection happens at compile
  time via `#ifdef`. There is no runtime CPUID dispatch; you must compile with
  the appropriate flags (`-mcpu=native` for ARM, `-mavx2 -mfma` for x86 AVX2,
  `-mavx512f` for AVX512) to activate the relevant path. A single binary cannot
  auto-select.
- **IV solver convergence is not guaranteed for extreme inputs.** Very deep OTM
  options with near-zero market prices can have vanishing vega, causing the
  solver to exit early with a best-effort estimate rather than a converged result.
- **Windows is untested.** `_Alignas` requires C11; `aligned_alloc` requires a
  POSIX/glibc environment. The benchmark and example code uses POSIX-only APIs.

## Quick Start

### Prerequisites

- C compiler: Clang or GCC
- CMake 3.16+
- A supported platform: ARM64, x86-64, or any C99 scalar target

### Installation

```bash
git clone https://github.com/NotoriousLOB/wondervol.git
cd wondervol
mkdir build && cd build
cmake ..
make -j$(nproc)
sudo make install
```

### Basic Usage

```c
#include <wondervol.h>
#include <stdio.h>

int main(void) {
    // Price a European call option
    double price = bsm_price(
        100.0,  // spot
        105.0,  // strike
        0.5,    // time to expiry (years)
        0.05,   // risk-free rate
        0.02,   // dividend yield
        0.25,   // volatility
        1       // 1=call, 0=put
    );
    printf("Call price: $%.4f\n", price);

    // Recover implied volatility from a market price
    double iv = bsm_implied_vol(100.0, 105.0, 0.5, 0.05, 0.02, 8.50, 1);
    printf("Implied volatility: %.2f%%\n", iv * 100.0);

    return 0;
}
```

Compile:

```bash
# ARM64
clang -O3 -mcpu=native -I/usr/local/include my_program.c -o my_program -lm

# x86-64 with AVX2
clang -O3 -mavx2 -mfma -I/usr/local/include my_program.c -o my_program -lm

# x86-64 with AVX512
clang -O3 -mavx512f -mavx512dq -I/usr/local/include my_program.c -o my_program -lm

# Any platform (scalar fallback)
clang -O3 -I/usr/local/include my_program.c -o my_program -lm
```

## Building

### CMake Options

| Option | Default | Description |
|--------|---------|-------------|
| `WV_BUILD_TESTS` | ON | Build test suite |
| `WV_BUILD_BENCHMARKS` | ON | Build benchmarks |
| `WV_BUILD_EXAMPLES` | ON | Build examples |
| `WV_ENABLE_OPENMP` | ON | Enable OpenMP parallelization |
| `WV_STRICT_TRAP` | OFF | Abort on validation failures |

```bash
# Full build
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)

# Library only
cmake .. -DWV_BUILD_TESTS=OFF -DWV_BUILD_BENCHMARKS=OFF -DWV_BUILD_EXAMPLES=OFF
make

# Debug build with strict validation
cmake .. -DCMAKE_BUILD_TYPE=Debug -DWV_STRICT_TRAP=ON
make
```

## Running Tests

```bash
cd build
ctest --output-on-failure

# Individual binaries
./tests/test_bsm          # BSM pricing and CDF accuracy
./tests/test_iv           # IV solver, all vols, deep OTM
./tests/test_dispersion   # Implied correlation and scanner
./tests/test_platform     # SIMD backend detection, vec vs scalar agreement
```

## Benchmarks

```bash
cd build

# IV throughput (100k options)
./benchmarks/bench_iv

# SIMD vs scalar comparison
./benchmarks/bench_simd_scalar 100000 10

# Dispersion scanner
./benchmarks/bench_dispersion 1000 50
```

### Example Output (ARM64, NEON)

```
=== Results ===
SIMD backend  : NEON (2-wide)
Scalar (ST)   : 28.9 ms  (3,460,000 opts/sec)
Scalar (OMP)  : 12.1 ms  (8,264,000 opts/sec)
SIMD + OMP    : 12.3 ms  (8,100,000 opts/sec)

Speedup (SIMD vs Scalar ST)  : 2.35x
Speedup (SIMD vs Scalar OMP) : 0.98x
Speedup (OMP vs ST)          : 2.39x
```

The NEON 2-wide path's speedup over scalar is modest because `exp`/`log` are evaluated per-lane in scalar (see Limitations). The OpenMP parallelism provides the majority of the throughput gain on multi-core hardware.

## Dispersion Trading Example

```c
#include <wondervol.h>
#include <stdio.h>

int main(void) {
    ComponentVol basket[] = {
        {0.25, 0.30, {0}},  // 25% weight, 30% vol
        {0.20, 0.28, {0}},
        {0.15, 0.35, {0}},
        {0.12, 0.32, {0}},
        {0.10, 0.25, {0}},
        {0.08, 0.27, {0}},
        {0.10, 0.29, {0}},
    };
    size_t n_comps = 7;

    double etf_iv = 0.22;
    double rho = implied_correlation(etf_iv, n_comps, basket);
    printf("Implied correlation: %.4f\n", rho);

    if (rho > 0.70)
        printf("Signal: SELL ETF vol (long dispersion)\n");
    else if (rho < 0.40)
        printf("Signal: BUY ETF vol (short dispersion)\n");

    return 0;
}
```

## API Reference

The arch-neutral names below work on all platforms. The `_arm` suffixed names are aliases on NEON platforms and remain available for backward compatibility.

### Pricing

```c
double bsm_price(double S, double K, double tau, double r, double q,
                 double sigma, int is_call);
```

Returns the BSM European option price. Returns intrinsic value when `tau <= 0`.

### Implied Volatility (scalar)

```c
double bsm_implied_vol(double S, double K, double tau, double r, double q,
                       double mkt_price, int is_call);
```

Returns IV in `[0.001, 5.0]`. Returns `0.0` when `tau <= 0`. Never crashes on bad inputs — returns a clamped best-effort result.

### Implied Volatility (vectorized array)

```c
void bsm_implied_vol_vec(const double* S, const double* K, const double* tau,
                         const double* r, const double* q,
                         const double* mkt_price, const uint8_t* is_call,
                         size_t n, double* iv_out);
```

Dispatches to the widest available SIMD path. All pointers must be non-null and non-overlapping (`restrict`). `n` may be any size; remainders fall to scalar.

Platform-specific names: `bsm_implied_vol_vec_arm`, `bsm_implied_vol_vec_x86`, `bsm_implied_vol_vec_scalar`.

### Implied Correlation

```c
double implied_correlation(double etf_iv, size_t n_comps,
                           const ComponentVol* comps);
```

Returns `ρ ∈ [−1, 1]`, `NAN` for empty baskets, `1.0` for single-component baskets.

### Dispersion Scanner

```c
void scan_dispersion_opps(const double* etf_ivs, const double* etf_vegas,
                          size_t n_etf, const ComponentVol* comps,
                          size_t n_comps, double target_rho, double thresh,
                          double basket_vega, DispersionOpportunity* opps);
```

Writes into `opps[n_etf]`. Entries where `score <= thresh` have `score = 0.0` and other fields are undefined.

## Data Structures

```c
typedef struct {
    double weight;             // Component weight in basket (should sum to 1)
    double iv;                 // Implied volatility
    _Alignas(64) double pad[6]; // Cache-line padding for SIMD prefetch
} ComponentVol;

typedef struct {
    double implied_rho;             // Implied correlation extracted from ETF IV
    double score;                   // |rho − target_rho|; 0 if below threshold
    double etf_vega;                // ETF vega at current market
    double suggested_etf_notional;  // −basket_vega / etf_vega (hedge ratio)
    double basket_vega;             // Input basket vega (passed through)
} DispersionOpportunity;
```

## Performance

Benchmarks on ARM64 (1.5GHz potato, NEON, no OpenMP available on this build):

| Operation | Backend | Throughput |
|-----------|---------|-----------|
| IV solve (100k options) | ARM NEON | ~8M opts/sec |
| IV solve (100k options) | Scalar ST | ~3.5M opts/sec |
| Dispersion scan (1k ETFs, 50 comps) | Scalar | <1 ms |

## FAQ

**Q:** Why no Heston / rough vol / local vol / whatever stochastic circle-jerk is trending on Twitter this week?

**A:** Because 99% of the time you don't need it. We love academic masturbation that runs at 12 Hz as much as the next guy, but most of the code we ship is for losing money _fast_.


**Q:** Why header-only?

**A:** Because linking is for people who enjoy waiting on the linker while their code is already 3x slower than it needs to be.


**Q:** Can I use this in production?

**A:** Absolutely! But only if you have the desire to see your book torn asunder at the tentacles of an angery gamma kraken. We'd appreciate a heads-up so we can watch.


**Q:** My boss wants Python bindings.

**A:** Tell your boss to learn C or keep paying the numpy tax. We don’t do that here.

## Acknowledgements

- Black, F. & Scholes, M. (1973). *The pricing of options and corporate
  liabilities.* Journal of Political Economy 81(3):637–654.
- Merton, R.C. (1973). *Theory of rational option pricing.* Bell Journal of
  Economics 4(1):141–183.
- Corrado, C.J. & Miller, T.W. (1996). *A note on a simple, accurate formula
  to compute implied standard deviations.* Journal of Banking & Finance
  20(3):595–603. *(IV initial guess)*
- Sun Microsystems / FreeBSD libm `e_erf.c` — FDLIBM piecewise rational erfc
  used for the normal CDF. *(~1 ulp accuracy)*
- ARM NEON, Intel AVX2/AVX512 intrinsics documentation.

