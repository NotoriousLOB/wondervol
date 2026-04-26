#pragma once
#ifndef WONDERVOL_H
#define WONDERVOL_H

#include <math.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <stdint.h>
#include <errno.h>
#ifdef _OPENMP
#   include <omp.h>
#endif

/* ===================================================================
 * Platform SIMD detection
 * =================================================================== */
#if defined(__ARM_NEON) || defined(__ARM_NEON__)
#   define WV_ARCH_NEON 1
#   include <arm_neon.h>
#elif defined(__AVX512F__)
#   define WV_ARCH_AVX512 1
#   include <immintrin.h>
#elif defined(__AVX2__)
#   define WV_ARCH_AVX2 1
#   include <immintrin.h>
#else
#   define WV_ARCH_SCALAR 1
#endif

/* ===================================================================
 * Validation macro
 * =================================================================== */
#ifndef WV_STRICT_TRAP
#   define WV_STRICT_TRAP 0
#endif

#define WV_REQUIRE(cond, msg) do { \
    if (!(cond)) { \
        static int warned = 0; \
        if (!warned) { \
            fprintf(stderr, "WONDERVOL_VIOLATION: %s (%s:%d) — %s\n", \
                    #cond, __FILE__, __LINE__, msg); \
            warned = 1; \
        } \
        if (WV_STRICT_TRAP) abort(); \
        errno = EINVAL; \
    } \
} while (0)

#ifndef WV_OMP_CHUNK_SIZE
#define WV_OMP_CHUNK_SIZE 256
#endif

/* ===================================================================
 * Constants
 * =================================================================== */
static const double WV_SQRT2        = 1.41421356237309504880168872420969808;
static const double WV_INV_SQRT_2PI = 0.39894228040143267793994605993438186;

/*
 * FDLIBM precise_erfc — piecewise rational approximation, ~1 ulp accuracy.
 * Ported from erf.h (Sun Microsystems / FreeBSD libm heritage).
 * Four regions: |x|<0.84375, |x|<1.25, |x|<26.5, |x|>=26.5.
 */
static const double wv_erx  = 8.45062911510467529297e-01;
static const double wv_pp0  =  1.28379167095512558561e-01;
static const double wv_pp1  = -3.25042107247001465170e-01;
static const double wv_pp2  = -2.84817495755985104766e-02;
static const double wv_pp3  = -5.77027029648944159157e-03;
static const double wv_pp4  = -2.37630166566501626084e-05;
static const double wv_qq1  =  3.97917223959155352819e-01;
static const double wv_qq2  =  6.50222499887672944485e-02;
static const double wv_qq3  =  5.08130628187576562776e-03;
static const double wv_qq4  =  1.32494738004321644526e-04;
static const double wv_qq5  = -3.96022827877536823420e-06;
static const double wv_pa0  = -2.36211856075265944077e-03;
static const double wv_pa1  =  4.14856118683748331466e-01;
static const double wv_pa2  = -3.72207876035701323847e-01;
static const double wv_pa3  =  3.18346619901161753674e-01;
static const double wv_pa4  = -1.10894694282396677476e-01;
static const double wv_pa5  =  3.54783043256182359371e-02;
static const double wv_pa6  = -2.16637559486879084300e-03;
static const double wv_qa1  =  1.06420880400844228286e-01;
static const double wv_qa2  =  5.40397917702171048937e-01;
static const double wv_qa3  =  7.18286544141962662268e-02;
static const double wv_qa4  =  1.26171219808761642112e-01;
static const double wv_qa5  =  1.36370839120290507362e-02;
static const double wv_qa6  =  1.19844998467991074170e-02;
static const double wv_ra0  = -9.86494403484714822705e-03;
static const double wv_ra1  = -6.93858572707181764372e-01;
static const double wv_ra2  = -1.05586262253232909814e+01;
static const double wv_ra3  = -6.23753324503260060396e+01;
static const double wv_ra4  = -1.62396669462573470355e+02;
static const double wv_ra5  = -1.84605092906711035994e+02;
static const double wv_ra6  = -8.12874355063065934246e+01;
static const double wv_ra7  = -9.81432934416914548592e+00;
static const double wv_sa1  =  1.96512716674392571292e+01;
static const double wv_sa2  =  1.37657754143519042600e+02;
static const double wv_sa3  =  4.34565877475292828876e+02;
static const double wv_sa4  =  6.45387271733244479411e+02;
static const double wv_sa5  =  4.29008140027567822286e+02;
static const double wv_sa6  =  1.08635005541778435161e+02;
static const double wv_sa7  =  6.57024977031928170135e+00;
static const double wv_sa8  = -6.04244152148580987438e-02;
static const double wv_rb0  = -9.86494292470009928597e-03;
static const double wv_rb1  = -7.99283237680523066574e-01;
static const double wv_rb2  = -1.77579549177547519889e+01;
static const double wv_rb3  = -1.60636384855821916062e+02;
static const double wv_rb4  = -6.37566443368389686328e+02;
static const double wv_rb5  = -1.02509513161107724954e+03;
static const double wv_rb6  = -4.83519191608651397019e+02;
static const double wv_sb1  =  3.03380607434824582924e+01;
static const double wv_sb2  =  3.25792512996573918826e+02;
static const double wv_sb3  =  1.53672958608443695914e+03;
static const double wv_sb4  =  3.19985821950859553918e+03;
static const double wv_sb5  =  2.55305040643316442583e+03;
static const double wv_sb6  =  4.74528541206955367215e+02;
static const double wv_sb7  = -2.24409524465858183362e+01;

static inline double wv_precise_erfc(double x) {
    if (x != x) return x + x;
    if (x == 0.0) return 1.0;
    const double ax = fabs(x);
    if (ax < 0.84375) {
        if (ax < 0x1.0p-54) return 1.0 - x;
        const double z = ax * ax;
        double r = wv_pp4; r = r*z+wv_pp3; r = r*z+wv_pp2; r = r*z+wv_pp1; r = r*z+wv_pp0;
        double s = wv_qq5; s = s*z+wv_qq4; s = s*z+wv_qq3; s = s*z+wv_qq2; s = s*z+wv_qq1; s = s*z+1.0;
        const double y = r / s;
        return 1.0 - x - x * y;
    }
    if (ax < 1.25) {
        const double t = ax - 1.0;
        double P = wv_pa6; P=P*t+wv_pa5; P=P*t+wv_pa4; P=P*t+wv_pa3; P=P*t+wv_pa2; P=P*t+wv_pa1; P=P*t+wv_pa0;
        double Q = wv_qa6; Q=Q*t+wv_qa5; Q=Q*t+wv_qa4; Q=Q*t+wv_qa3; Q=Q*t+wv_qa2; Q=Q*t+wv_qa1; Q=Q*t+1.0;
        const double erfc_ax = 1.0 - wv_erx - P / Q;
        return (x >= 0.0) ? erfc_ax : 2.0 - erfc_ax;
    }
    if (ax < 26.5) {
        const double s = 1.0 / (ax * ax);
        double R, S;
        if (ax < 2.85711669921875) {
            R = wv_ra7; R=R*s+wv_ra6; R=R*s+wv_ra5; R=R*s+wv_ra4; R=R*s+wv_ra3; R=R*s+wv_ra2; R=R*s+wv_ra1; R=R*s+wv_ra0;
            S = wv_sa8; S=S*s+wv_sa7; S=S*s+wv_sa6; S=S*s+wv_sa5; S=S*s+wv_sa4; S=S*s+wv_sa3; S=S*s+wv_sa2; S=S*s+wv_sa1; S=S*s+1.0;
        } else {
            R = wv_rb6; R=R*s+wv_rb5; R=R*s+wv_rb4; R=R*s+wv_rb3; R=R*s+wv_rb2; R=R*s+wv_rb1; R=R*s+wv_rb0;
            S = wv_sb7; S=S*s+wv_sb6; S=S*s+wv_sb5; S=S*s+wv_sb4; S=S*s+wv_sb3; S=S*s+wv_sb2; S=S*s+wv_sb1; S=S*s+1.0;
        }
        const double tail = exp(-ax*ax - 0.5625 + R/S) / ax;
        return (x >= 0.0) ? tail : 2.0 - tail;
    }
    return (x >= 0.0) ? 0.0 : 2.0;
}

/* ===================================================================
 * Scalar CDF / PDF — Φ(x) = erfc(-x/√2) / 2
 * =================================================================== */
static inline double norm_pdf(double x) {
    return WV_INV_SQRT_2PI * exp(-0.5 * x * x);
}

static inline double norm_cdf(double x) {
    return 0.5 * wv_precise_erfc(-x / WV_SQRT2);
}

/* ===================================================================
 * Corrado-Miller (1996) IV initial guess
 * Reference: Corrado & Miller, J. Banking & Finance 20:151–169.
 *
 * σ₀ = (1/Se^{-qτ}) × √(2π/τ) × [C̃ - (F-K)/2 + √((C̃-(F-K)/2)² - (F-K)²/π)]
 * where F = Se^{(r-q)τ}, C̃ is the undiscounted call price.
 * =================================================================== */
static inline double wv_iv_initial_guess(double S, double K, double tau, double r,
                                         double q, double mkt_price, int is_call) {
    double disc_r = exp(-r * tau);
    double F      = S * exp((r - q) * tau);
    double S_q    = S * exp(-q * tau);

    /* Convert put to undiscounted call via parity: C̃ = P/disc_r + (F-K) */
    double C_undisc = is_call ? (mkt_price / disc_r) : (mkt_price / disc_r + (F - K));

    double mid   = C_undisc - 0.5 * (F - K);
    double disc2 = (F - K) * (F - K) / M_PI;
    double inner = mid * mid - disc2;
    if (inner < 0.0) inner = 0.0;

    double sigma = (sqrt(2.0 * M_PI / tau) / S_q) * (mid + sqrt(inner));
    if (!isfinite(sigma) || sigma < 0.01) sigma = 0.2;
    return fmax(0.001, fmin(5.0, sigma));
}

/* ===================================================================
 * Scalar BSM price
 * =================================================================== */
static inline double bsm_price(double S, double K, double tau, double r, double q,
                               double sigma, int is_call) {
    if (tau <= 0.0 || !isfinite(tau)) {
        return is_call ? fmax(S - K, 0.0) : fmax(K - S, 0.0);
    }
    double sqrt_tau = sqrt(tau);
    double d1 = (log(S/K) + (r - q + 0.5*sigma*sigma)*tau) / (sigma * sqrt_tau);
    double d2 = d1 - sigma * sqrt_tau;
    double disc_r = exp(-r * tau);
    double disc_q = exp(-q * tau);
    double Nd1 = norm_cdf(d1);
    double Nd2 = norm_cdf(d2);
    double call = S * disc_q * Nd1 - K * disc_r * Nd2;
    double put  = K * disc_r * (1.0 - Nd2) - S * disc_q * (1.0 - Nd1);
    return is_call ? call : put;
}

/* ===================================================================
 * Scalar IV solver — Newton-Raphson with Corrado-Miller initial guess
 * Convergence: |price diff| < WV_NR_TOL = 1e-12 or vega < WV_NR_TOL
 * =================================================================== */
#define WV_NR_TOL 1e-12

static inline double bsm_implied_vol(double S, double K, double tau, double r, double q,
                                     double mkt_price, int is_call) {
    WV_REQUIRE(S > 0 && K > 0 && isfinite(mkt_price), "invalid BSM params");
    if (tau <= 0.0) return 0.0;

    double sigma = wv_iv_initial_guess(S, K, tau, r, q, mkt_price, is_call);
    double sqrt_tau = sqrt(tau);
    double disc_q   = exp(-q * tau);
    double fwd_S    = S * disc_q;

    for (int i = 0; i < 100; ++i) {
        double price = bsm_price(S, K, tau, r, q, sigma, is_call);
        double diff  = price - mkt_price;
        if (fabs(diff) < WV_NR_TOL) return sigma;
        double d1   = (log(S/K) + (r - q + 0.5*sigma*sigma)*tau) / (sigma * sqrt_tau);
        double vega = fwd_S * norm_pdf(d1) * sqrt_tau;
        if (vega < WV_NR_TOL) break;
        double step = diff / vega;
        step  = fmax(fmin(step, 0.5 * sigma), -0.5 * sigma);
        sigma -= step;
        sigma = fmax(0.001, fmin(5.0, sigma));
    }
    return sigma;
}

/* ===================================================================
 * Data structures
 * =================================================================== */
typedef struct {
    double weight;
    double iv;
    _Alignas(64) double pad[6];
} ComponentVol;

typedef struct {
    double implied_rho;
    double score;
    double etf_vega;
    double suggested_etf_notional;
    double basket_vega;
} DispersionOpportunity;

/* ===================================================================
 * Scalar fallback array solvers — used when no SIMD is available,
 * also serve as arch-neutral reference implementations.
 * =================================================================== */
static inline void bsm_implied_vol_vec_scalar(const double* restrict S,
                                              const double* restrict K,
                                              const double* restrict tau,
                                              const double* restrict r,
                                              const double* restrict q,
                                              const double* restrict mkt_price,
                                              const uint8_t* restrict is_call,
                                              size_t n,
                                              double* restrict iv_out) {
    WV_REQUIRE(S && K && tau && r && q && mkt_price && is_call && iv_out, "non-null pointers");
    for (size_t i = 0; i < n; ++i) {
        iv_out[i] = bsm_implied_vol(S[i], K[i], tau[i], r[i], q[i], mkt_price[i], is_call[i]);
    }
}

static inline double implied_correlation_scalar(double etf_iv, size_t n_comps,
                                                const ComponentVol* restrict comps) {
    WV_REQUIRE(etf_iv > 0.0 && comps, "valid args");
    if (n_comps == 0) return NAN;
    if (n_comps == 1) return 1.0;
    double sum_w2_sig2 = 0.0, sum_w_sig = 0.0;
    for (size_t i = 0; i < n_comps; ++i) {
        double w = comps[i].weight, sig = comps[i].iv;
        sum_w2_sig2 += w * w * sig * sig;
        sum_w_sig   += w * sig;
    }
    double denom = sum_w_sig * sum_w_sig - sum_w2_sig2;
    if (denom < 1e-12) return 1.0;
    double rho = (etf_iv * etf_iv - sum_w2_sig2) / denom;
    return fmax(-1.0, fmin(1.0, rho));
}

static inline void scan_dispersion_opps_scalar(
    const double* restrict etf_ivs, const double* restrict etf_vegas, size_t n_etf,
    const ComponentVol* restrict comps, size_t n_comps,
    double target_rho, double thresh, double basket_vega,
    DispersionOpportunity* restrict opps) {

    WV_REQUIRE(etf_ivs && etf_vegas && comps && opps && n_etf > 0, "valid args");

    if (n_comps == 1) {
        for (size_t i = 0; i < n_etf; ++i) {
            opps[i].implied_rho = 1.0;
            opps[i].score = fabs(1.0 - target_rho);
            opps[i].etf_vega = etf_vegas[i];
            opps[i].suggested_etf_notional = (etf_vegas[i] > 1e-12) ? (-basket_vega / etf_vegas[i]) : 0.0;
            opps[i].basket_vega = basket_vega;
        }
        return;
    }

    double sum_w2_sig2 = 0.0, sum_w_sig = 0.0;
    for (size_t i = 0; i < n_comps; ++i) {
        double w = comps[i].weight, sig = comps[i].iv;
        sum_w2_sig2 += w * w * sig * sig;
        sum_w_sig   += w * sig;
    }
    double denom = sum_w_sig * sum_w_sig - sum_w2_sig2;
    if (denom < 1e-12) denom = 1.0;

    for (size_t i = 0; i < n_etf; ++i) {
        double rho = (etf_ivs[i] * etf_ivs[i] - sum_w2_sig2) / denom;
        double score = fabs(rho - target_rho);
        double vega_etf = etf_vegas[i];
        double notional = (vega_etf > 1e-12) ? (-basket_vega / vega_etf) : 0.0;
        if (score > thresh) {
            opps[i].implied_rho            = rho;
            opps[i].score                  = score;
            opps[i].etf_vega               = vega_etf;
            opps[i].suggested_etf_notional = notional;
            opps[i].basket_vega            = basket_vega;
        } else {
            opps[i].score = 0.0;
        }
    }
}

/* ===================================================================
 * ARM NEON backend (2-wide float64x2_t)
 * =================================================================== */
#if defined(WV_ARCH_NEON)

static inline float64x2_t wv_neon_exp_f64(float64x2_t x) {
    double v0 = exp(vgetq_lane_f64(x, 0));
    double v1 = exp(vgetq_lane_f64(x, 1));
    return vsetq_lane_f64(v1, vdupq_n_f64(v0), 1);
}

static inline float64x2_t wv_neon_log_f64(float64x2_t x) {
    double v0 = log(vgetq_lane_f64(x, 0));
    double v1 = log(vgetq_lane_f64(x, 1));
    return vsetq_lane_f64(v1, vdupq_n_f64(v0), 1);
}

static inline float64x2_t norm_pdf_neon(float64x2_t x) {
    float64x2_t xx = vmulq_f64(x, x);
    return vmulq_f64(vdupq_n_f64(WV_INV_SQRT_2PI),
                     wv_neon_exp_f64(vmulq_f64(vdupq_n_f64(-0.5), xx)));
}

/* norm_cdf_neon: per-lane scalar call into wv_precise_erfc.
 * NEON has no native fp64 erfc; lane extraction + scalar is correct and fast
 * enough since CDF is called O(50) times per IV solve, not O(n) per option.
 */
static inline float64x2_t norm_cdf_neon(float64x2_t x) {
    double c0 = norm_cdf(vgetq_lane_f64(x, 0));
    double c1 = norm_cdf(vgetq_lane_f64(x, 1));
    return vsetq_lane_f64(c1, vdupq_n_f64(c0), 1);
}

static inline float64x2_t bsm_price_neon(float64x2_t S, float64x2_t K, float64x2_t tau,
                                         float64x2_t r, float64x2_t q, float64x2_t sigma,
                                         uint64x2_t call_mask) {
    float64x2_t zero = vdupq_n_f64(0.0);
    float64x2_t one  = vdupq_n_f64(1.0);

    float64x2_t call_int = vmaxq_f64(vsubq_f64(S, K), zero);
    float64x2_t put_int  = vmaxq_f64(vsubq_f64(K, S), zero);
    float64x2_t expired_price = vbslq_f64(call_mask, call_int, put_int);

    uint64x2_t expired = vcleq_f64(tau, zero);

    if (vgetq_lane_u64(expired, 0) && vgetq_lane_u64(expired, 1)) {
        return expired_price;
    }

    float64x2_t safe_tau    = vmaxq_f64(tau, vdupq_n_f64(1e-30));
    float64x2_t sqrt_tau    = vsqrtq_f64(safe_tau);
    float64x2_t sig_sqrt_tau = vmulq_f64(sigma, sqrt_tau);
    float64x2_t logSK       = wv_neon_log_f64(vdivq_f64(S, K));
    float64x2_t half_sig2   = vmulq_f64(vdupq_n_f64(0.5), vmulq_f64(sigma, sigma));
    float64x2_t drift       = vfmaq_f64(logSK, vsubq_f64(vaddq_f64(r, half_sig2), q), safe_tau);
    float64x2_t d1          = vdivq_f64(drift, sig_sqrt_tau);
    float64x2_t d2          = vsubq_f64(d1, sig_sqrt_tau);

    float64x2_t Nd1   = norm_cdf_neon(d1);
    float64x2_t Nd2   = norm_cdf_neon(d2);
    float64x2_t disc_r = wv_neon_exp_f64(vnegq_f64(vmulq_f64(r, safe_tau)));
    float64x2_t disc_q = wv_neon_exp_f64(vnegq_f64(vmulq_f64(q, safe_tau)));
    float64x2_t fwd_S  = vmulq_f64(S, disc_q);
    float64x2_t fwd_K  = vmulq_f64(K, disc_r);

    float64x2_t call_pr = vsubq_f64(vmulq_f64(fwd_S, Nd1), vmulq_f64(fwd_K, Nd2));
    float64x2_t put_pr  = vsubq_f64(vmulq_f64(fwd_K, vsubq_f64(one, Nd2)),
                                     vmulq_f64(fwd_S, vsubq_f64(one, Nd1)));
    float64x2_t live_price = vbslq_f64(call_mask, call_pr, put_pr);

    return vbslq_f64(expired, expired_price, live_price);
}

static inline void bsm_implied_vol_vec_arm(const double* restrict S,
                                           const double* restrict K,
                                           const double* restrict tau,
                                           const double* restrict r,
                                           const double* restrict q,
                                           const double* restrict mkt_price,
                                           const uint8_t* restrict is_call,
                                           size_t n,
                                           double* restrict iv_out) {
    WV_REQUIRE(S && K && tau && r && q && mkt_price && is_call && iv_out, "non-null pointers");

    #ifdef _OPENMP
    #pragma omp parallel for schedule(static, WV_OMP_CHUNK_SIZE)
    #endif
    for (size_t i = 0; i < (n & ~1ULL); i += 2) {
        float64x2_t Sv   = vld1q_f64(S + i);
        float64x2_t Kv   = vld1q_f64(K + i);
        float64x2_t tauv = vld1q_f64(tau + i);
        float64x2_t rv   = vld1q_f64(r + i);
        float64x2_t qv   = vld1q_f64(q + i);
        float64x2_t mv   = vld1q_f64(mkt_price + i);

        uint64x2_t call_mask = vdupq_n_u64(0);
        call_mask = vsetq_lane_u64(is_call[i]   ? ~0ULL : 0ULL, call_mask, 0);
        call_mask = vsetq_lane_u64(is_call[i+1] ? ~0ULL : 0ULL, call_mask, 1);

        float64x2_t one  = vdupq_n_f64(1.0);

        /* Corrado-Miller initial guess, per lane */
        double sg0 = wv_iv_initial_guess(S[i],   K[i],   tau[i],   r[i],   q[i],   mkt_price[i],   is_call[i]);
        double sg1 = wv_iv_initial_guess(S[i+1], K[i+1], tau[i+1], r[i+1], q[i+1], mkt_price[i+1], is_call[i+1]);
        float64x2_t sigma = vsetq_lane_f64(sg1, vdupq_n_f64(sg0), 1);

        float64x2_t safe_tau = vmaxq_f64(tauv, vdupq_n_f64(1e-30));
        float64x2_t sqrt_tau = vsqrtq_f64(safe_tau);
        float64x2_t logSK    = wv_neon_log_f64(vdivq_f64(Sv, Kv));
        float64x2_t disc_r   = wv_neon_exp_f64(vnegq_f64(vmulq_f64(rv, safe_tau)));
        float64x2_t disc_q   = wv_neon_exp_f64(vnegq_f64(vmulq_f64(qv, safe_tau)));
        float64x2_t fwd_S    = vmulq_f64(Sv, disc_q);
        float64x2_t fwd_K    = vmulq_f64(Kv, disc_r);

        /* Lane-disable mask: lanes where vega < WV_NR_TOL are frozen */
        uint64x2_t active = vdupq_n_u64(~0ULL);

        for (int iter = 0; iter < 50; ++iter) {
            float64x2_t halfv2 = vmulq_f64(vdupq_n_f64(0.5), vmulq_f64(sigma, sigma));
            float64x2_t num    = vfmaq_f64(logSK, vsubq_f64(vaddq_f64(rv, halfv2), qv), safe_tau);
            float64x2_t d1     = vdivq_f64(num, vmulq_f64(sigma, sqrt_tau));
            float64x2_t d2     = vsubq_f64(d1, vmulq_f64(sigma, sqrt_tau));

            float64x2_t Nd1   = norm_cdf_neon(d1);
            float64x2_t Nd2   = norm_cdf_neon(d2);
            float64x2_t call_pr = vsubq_f64(vmulq_f64(fwd_S, Nd1), vmulq_f64(fwd_K, Nd2));
            float64x2_t put_pr  = vsubq_f64(vmulq_f64(fwd_K, vsubq_f64(one, Nd2)),
                                            vmulq_f64(fwd_S, vsubq_f64(one, Nd1)));
            float64x2_t price  = vbslq_f64(call_mask, call_pr, put_pr);
            float64x2_t diff   = vsubq_f64(price, mv);
            float64x2_t abs_diff = vabsq_f64(diff);

            if (vgetq_lane_f64(abs_diff, 0) < WV_NR_TOL && vgetq_lane_f64(abs_diff, 1) < WV_NR_TOL) {
                break;
            }

            float64x2_t vega = vmulq_f64(fwd_S, vmulq_f64(norm_pdf_neon(d1), sqrt_tau));

            /* Freeze lanes where vega is singular */
            uint64x2_t vega_ok = vcgeq_f64(vega, vdupq_n_f64(WV_NR_TOL));
            active = vandq_u64(active, vega_ok);
            if (!vgetq_lane_u64(active, 0) && !vgetq_lane_u64(active, 1)) break;

            float64x2_t step       = vdivq_f64(diff, vmaxq_f64(vega, vdupq_n_f64(WV_NR_TOL)));
            float64x2_t half_sigma = vmulq_f64(vdupq_n_f64(0.5), sigma);
            step = vminq_f64(step, half_sigma);
            step = vmaxq_f64(step, vnegq_f64(half_sigma));

            /* Only update active lanes */
            float64x2_t new_sigma = vsubq_f64(sigma, step);
            new_sigma = vmaxq_f64(vdupq_n_f64(0.001), vminq_f64(vdupq_n_f64(5.0), new_sigma));
            sigma = vbslq_f64(active, new_sigma, sigma);
        }

        vst1q_f64(iv_out + i, sigma);

        if (vgetq_lane_f64(tauv, 0) <= 0.0) iv_out[i]     = 0.0;
        if (vgetq_lane_f64(tauv, 1) <= 0.0) iv_out[i + 1] = 0.0;
    }

    /* Scalar tail */
    if (n % 2 != 0) {
        size_t last = n - 1;
        iv_out[last] = bsm_implied_vol(S[last], K[last], tau[last], r[last], q[last],
                                       mkt_price[last], is_call[last]);
    }
}

static inline double implied_correlation_arm(double etf_iv, size_t n_comps,
                                             const ComponentVol* restrict comps) {
    return implied_correlation_scalar(etf_iv, n_comps, comps);
}

static inline void scan_dispersion_opps_arm(
    const double* restrict etf_ivs, const double* restrict etf_vegas, size_t n_etf,
    const ComponentVol* restrict comps, size_t n_comps,
    double target_rho, double thresh, double basket_vega,
    DispersionOpportunity* restrict opps) {

    WV_REQUIRE(etf_ivs && etf_vegas && comps && opps && n_etf > 0, "valid args");

    if (n_comps == 1) {
        #ifdef _OPENMP
        #pragma omp parallel for schedule(static)
        #endif
        for (size_t i = 0; i < n_etf; ++i) {
            opps[i].implied_rho = 1.0;
            opps[i].score = fabs(1.0 - target_rho);
            opps[i].etf_vega = etf_vegas[i];
            opps[i].suggested_etf_notional = (etf_vegas[i] > 1e-12) ? (-basket_vega / etf_vegas[i]) : 0.0;
            opps[i].basket_vega = basket_vega;
        }
        return;
    }

    double sum_w2_sig2 = 0.0, sum_w_sig = 0.0;
    for (size_t i = 0; i < n_comps; ++i) {
        double w = comps[i].weight, sig = comps[i].iv;
        sum_w2_sig2 += w * w * sig * sig;
        sum_w_sig   += w * sig;
    }
    double denom = sum_w_sig * sum_w_sig - sum_w2_sig2;
    if (denom < 1e-12) denom = 1.0;

    #ifdef _OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for (size_t i = 0; i < n_etf; ++i) {
        double rho     = (etf_ivs[i] * etf_ivs[i] - sum_w2_sig2) / denom;
        double score   = fabs(rho - target_rho);
        double vega_etf = etf_vegas[i];
        double notional = (vega_etf > 1e-12) ? (-basket_vega / vega_etf) : 0.0;
        if (score > thresh) {
            opps[i].implied_rho            = rho;
            opps[i].score                  = score;
            opps[i].etf_vega               = vega_etf;
            opps[i].suggested_etf_notional = notional;
            opps[i].basket_vega            = basket_vega;
        } else {
            opps[i].score = 0.0;
        }
    }
}

/* ===================================================================
 * AVX2 backend (4-wide __m256d)
 * =================================================================== */
#elif defined(WV_ARCH_AVX2)

static inline __m256d wv_avx2_exp_f64(__m256d x) {
    double v[4];
    _mm256_storeu_pd(v, x);
    v[0] = exp(v[0]); v[1] = exp(v[1]); v[2] = exp(v[2]); v[3] = exp(v[3]);
    return _mm256_loadu_pd(v);
}

static inline __m256d wv_avx2_log_f64(__m256d x) {
    double v[4];
    _mm256_storeu_pd(v, x);
    v[0] = log(v[0]); v[1] = log(v[1]); v[2] = log(v[2]); v[3] = log(v[3]);
    return _mm256_loadu_pd(v);
}

static inline __m256d norm_pdf_avx2(__m256d x) {
    __m256d xx   = _mm256_mul_pd(x, x);
    __m256d neg_half = _mm256_set1_pd(-0.5);
    return _mm256_mul_pd(_mm256_set1_pd(WV_INV_SQRT_2PI),
                         wv_avx2_exp_f64(_mm256_mul_pd(neg_half, xx)));
}

/*
 * Cody Algorithm 715 erfc, AVX2 4-wide.
 * Uses scalar per-lane evaluation for correctness; region blending
 * is done via _mm256_blendv_pd.
 */
static inline __m256d norm_cdf_avx2(__m256d x) {
    double v[4], r[4];
    _mm256_storeu_pd(v, x);
    for (int i = 0; i < 4; i++) r[i] = norm_cdf(v[i]);
    return _mm256_loadu_pd(r);
}

static inline __m256d bsm_price_avx2(__m256d S, __m256d K, __m256d tau,
                                     __m256d r, __m256d q, __m256d sigma,
                                     __m256d call_mask) {
    __m256d zero     = _mm256_setzero_pd();
    __m256d one      = _mm256_set1_pd(1.0);
    __m256d call_int = _mm256_max_pd(_mm256_sub_pd(S, K), zero);
    __m256d put_int  = _mm256_max_pd(_mm256_sub_pd(K, S), zero);
    __m256d exp_pr   = _mm256_blendv_pd(put_int, call_int, call_mask);

    __m256d expired  = _mm256_cmp_pd(tau, zero, _CMP_LE_OQ);
    /* fast path: check if all 4 expired */
    if (_mm256_movemask_pd(expired) == 0xF) return exp_pr;

    __m256d safe_tau    = _mm256_max_pd(tau, _mm256_set1_pd(1e-30));
    __m256d sqrt_tau    = _mm256_sqrt_pd(safe_tau);
    __m256d sig_sqrt    = _mm256_mul_pd(sigma, sqrt_tau);
    __m256d logSK       = wv_avx2_log_f64(_mm256_div_pd(S, K));
    __m256d half_sig2   = _mm256_mul_pd(_mm256_set1_pd(0.5), _mm256_mul_pd(sigma, sigma));
    __m256d drift       = _mm256_fmadd_pd(_mm256_sub_pd(_mm256_add_pd(r, half_sig2), q),
                                          safe_tau, logSK);
    __m256d d1          = _mm256_div_pd(drift, sig_sqrt);
    __m256d d2          = _mm256_sub_pd(d1, sig_sqrt);

    __m256d Nd1   = norm_cdf_avx2(d1);
    __m256d Nd2   = norm_cdf_avx2(d2);
    __m256d disc_r = wv_avx2_exp_f64(_mm256_mul_pd(_mm256_set1_pd(-1.0), _mm256_mul_pd(r, safe_tau)));
    __m256d disc_q = wv_avx2_exp_f64(_mm256_mul_pd(_mm256_set1_pd(-1.0), _mm256_mul_pd(q, safe_tau)));
    __m256d fwd_S  = _mm256_mul_pd(S, disc_q);
    __m256d fwd_K  = _mm256_mul_pd(K, disc_r);

    __m256d call_pr = _mm256_sub_pd(_mm256_mul_pd(fwd_S, Nd1), _mm256_mul_pd(fwd_K, Nd2));
    __m256d put_pr  = _mm256_sub_pd(_mm256_mul_pd(fwd_K, _mm256_sub_pd(one, Nd2)),
                                    _mm256_mul_pd(fwd_S, _mm256_sub_pd(one, Nd1)));
    __m256d live    = _mm256_blendv_pd(put_pr, call_pr, call_mask);
    return _mm256_blendv_pd(live, exp_pr, expired);
}

static inline void bsm_implied_vol_vec_x86(const double* restrict S,
                                           const double* restrict K,
                                           const double* restrict tau,
                                           const double* restrict r,
                                           const double* restrict q,
                                           const double* restrict mkt_price,
                                           const uint8_t* restrict is_call,
                                           size_t n,
                                           double* restrict iv_out) {
    WV_REQUIRE(S && K && tau && r && q && mkt_price && is_call && iv_out, "non-null pointers");

    __m256d one = _mm256_set1_pd(1.0);

    #ifdef _OPENMP
    #pragma omp parallel for schedule(static, WV_OMP_CHUNK_SIZE)
    #endif
    for (size_t i = 0; i < (n & ~3ULL); i += 4) {
        __m256d Sv   = _mm256_loadu_pd(S + i);
        __m256d Kv   = _mm256_loadu_pd(K + i);
        __m256d tauv = _mm256_loadu_pd(tau + i);
        __m256d rv   = _mm256_loadu_pd(r + i);
        __m256d qv   = _mm256_loadu_pd(q + i);
        __m256d mv   = _mm256_loadu_pd(mkt_price + i);

        double cm[4] = {
            (double)(is_call[i]  ), (double)(is_call[i+1]),
            (double)(is_call[i+2]), (double)(is_call[i+3])
        };
        __m256d call_mask_f = _mm256_loadu_pd(cm);
        __m256d call_mask   = _mm256_cmp_pd(call_mask_f, _mm256_setzero_pd(), _CMP_NEQ_OQ);

        /* Corrado-Miller initial guess per lane */
        double sg[4];
        for (int j = 0; j < 4; j++)
            sg[j] = wv_iv_initial_guess(S[i+j], K[i+j], tau[i+j], r[i+j], q[i+j],
                                        mkt_price[i+j], is_call[i+j]);
        __m256d sigma = _mm256_loadu_pd(sg);

        __m256d safe_tau = _mm256_max_pd(tauv, _mm256_set1_pd(1e-30));
        __m256d sqrt_tau = _mm256_sqrt_pd(safe_tau);
        __m256d logSK    = wv_avx2_log_f64(_mm256_div_pd(Sv, Kv));
        __m256d disc_r   = wv_avx2_exp_f64(_mm256_mul_pd(_mm256_set1_pd(-1.0), _mm256_mul_pd(rv, safe_tau)));
        __m256d disc_q   = wv_avx2_exp_f64(_mm256_mul_pd(_mm256_set1_pd(-1.0), _mm256_mul_pd(qv, safe_tau)));
        __m256d fwd_S    = _mm256_mul_pd(Sv, disc_q);
        __m256d fwd_K    = _mm256_mul_pd(Kv, disc_r);

        __m256d active = _mm256_castsi256_pd(_mm256_set1_epi64x(-1));

        for (int iter = 0; iter < 50; ++iter) {
            __m256d halfv2 = _mm256_mul_pd(_mm256_set1_pd(0.5), _mm256_mul_pd(sigma, sigma));
            __m256d num    = _mm256_fmadd_pd(_mm256_sub_pd(_mm256_add_pd(rv, halfv2), qv),
                                             safe_tau, logSK);
            __m256d d1     = _mm256_div_pd(num, _mm256_mul_pd(sigma, sqrt_tau));
            __m256d d2     = _mm256_sub_pd(d1, _mm256_mul_pd(sigma, sqrt_tau));

            __m256d Nd1   = norm_cdf_avx2(d1);
            __m256d Nd2   = norm_cdf_avx2(d2);
            __m256d call_pr = _mm256_sub_pd(_mm256_mul_pd(fwd_S, Nd1), _mm256_mul_pd(fwd_K, Nd2));
            __m256d put_pr  = _mm256_sub_pd(_mm256_mul_pd(fwd_K, _mm256_sub_pd(one, Nd2)),
                                           _mm256_mul_pd(fwd_S, _mm256_sub_pd(one, Nd1)));
            __m256d price  = _mm256_blendv_pd(put_pr, call_pr, call_mask);
            __m256d diff   = _mm256_sub_pd(price, mv);
            __m256d abs_diff = _mm256_andnot_pd(_mm256_set1_pd(-0.0), diff);

            double ad[4];
            _mm256_storeu_pd(ad, abs_diff);
            if (ad[0] < WV_NR_TOL && ad[1] < WV_NR_TOL && ad[2] < WV_NR_TOL && ad[3] < WV_NR_TOL)
                break;

            __m256d vega    = _mm256_mul_pd(fwd_S, _mm256_mul_pd(norm_pdf_avx2(d1), sqrt_tau));
            __m256d vega_ok = _mm256_cmp_pd(vega, _mm256_set1_pd(WV_NR_TOL), _CMP_GE_OQ);
            active = _mm256_and_pd(active, vega_ok);
            if (_mm256_movemask_pd(active) == 0) break;

            __m256d step       = _mm256_div_pd(diff, _mm256_max_pd(vega, _mm256_set1_pd(WV_NR_TOL)));
            __m256d half_sigma = _mm256_mul_pd(_mm256_set1_pd(0.5), sigma);
            step = _mm256_min_pd(step, half_sigma);
            step = _mm256_max_pd(step, _mm256_sub_pd(_mm256_setzero_pd(), half_sigma));

            __m256d new_sigma = _mm256_sub_pd(sigma, step);
            new_sigma = _mm256_max_pd(_mm256_set1_pd(0.001), _mm256_min_pd(_mm256_set1_pd(5.0), new_sigma));
            sigma = _mm256_blendv_pd(sigma, new_sigma, active);
        }

        _mm256_storeu_pd(iv_out + i, sigma);

        double tv[4];
        _mm256_storeu_pd(tv, tauv);
        for (int j = 0; j < 4; j++)
            if (tv[j] <= 0.0) iv_out[i+j] = 0.0;
    }

    /* 4-wide remainder falls through to scalar */
    for (size_t i = (n & ~3ULL); i < n; i++) {
        iv_out[i] = bsm_implied_vol(S[i], K[i], tau[i], r[i], q[i], mkt_price[i], is_call[i]);
    }
}

static inline double implied_correlation_x86(double etf_iv, size_t n_comps,
                                             const ComponentVol* restrict comps) {
    return implied_correlation_scalar(etf_iv, n_comps, comps);
}

static inline void scan_dispersion_opps_x86(
    const double* restrict etf_ivs, const double* restrict etf_vegas, size_t n_etf,
    const ComponentVol* restrict comps, size_t n_comps,
    double target_rho, double thresh, double basket_vega,
    DispersionOpportunity* restrict opps) {
    scan_dispersion_opps_scalar(etf_ivs, etf_vegas, n_etf, comps, n_comps,
                                target_rho, thresh, basket_vega, opps);
}

/* ===================================================================
 * AVX512 backend (8-wide __m512d)
 * =================================================================== */
#elif defined(WV_ARCH_AVX512)

static inline __m512d wv_avx512_exp_f64(__m512d x) {
    double v[8];
    _mm512_storeu_pd(v, x);
    for (int i = 0; i < 8; i++) v[i] = exp(v[i]);
    return _mm512_loadu_pd(v);
}

static inline __m512d wv_avx512_log_f64(__m512d x) {
    double v[8];
    _mm512_storeu_pd(v, x);
    for (int i = 0; i < 8; i++) v[i] = log(v[i]);
    return _mm512_loadu_pd(v);
}

static inline __m512d norm_pdf_avx512(__m512d x) {
    __m512d xx = _mm512_mul_pd(x, x);
    return _mm512_mul_pd(_mm512_set1_pd(WV_INV_SQRT_2PI),
                         wv_avx512_exp_f64(_mm512_mul_pd(_mm512_set1_pd(-0.5), xx)));
}

static inline __m512d norm_cdf_avx512(__m512d x) {
    double v[8], r[8];
    _mm512_storeu_pd(v, x);
    for (int i = 0; i < 8; i++) r[i] = norm_cdf(v[i]);
    return _mm512_loadu_pd(r);
}

static inline __m512d bsm_price_avx512(__m512d S, __m512d K, __m512d tau,
                                       __m512d r, __m512d q, __m512d sigma,
                                       __mmask8 call_mask) {
    __m512d zero    = _mm512_setzero_pd();
    __m512d one     = _mm512_set1_pd(1.0);
    __m512d call_int = _mm512_max_pd(_mm512_sub_pd(S, K), zero);
    __m512d put_int  = _mm512_max_pd(_mm512_sub_pd(K, S), zero);
    __m512d exp_pr   = _mm512_mask_blend_pd(call_mask, put_int, call_int);

    __mmask8 expired = _mm512_cmp_pd_mask(tau, zero, _CMP_LE_OQ);
    if (expired == 0xFF) return exp_pr;

    __m512d safe_tau  = _mm512_max_pd(tau, _mm512_set1_pd(1e-30));
    __m512d sqrt_tau  = _mm512_sqrt_pd(safe_tau);
    __m512d sig_sqrt  = _mm512_mul_pd(sigma, sqrt_tau);
    __m512d logSK     = wv_avx512_log_f64(_mm512_div_pd(S, K));
    __m512d half_sig2 = _mm512_mul_pd(_mm512_set1_pd(0.5), _mm512_mul_pd(sigma, sigma));
    __m512d drift     = _mm512_fmadd_pd(_mm512_sub_pd(_mm512_add_pd(r, half_sig2), q),
                                        safe_tau, logSK);
    __m512d d1   = _mm512_div_pd(drift, sig_sqrt);
    __m512d d2   = _mm512_sub_pd(d1, sig_sqrt);

    __m512d Nd1   = norm_cdf_avx512(d1);
    __m512d Nd2   = norm_cdf_avx512(d2);
    __m512d disc_r = wv_avx512_exp_f64(_mm512_mul_pd(_mm512_set1_pd(-1.0), _mm512_mul_pd(r, safe_tau)));
    __m512d disc_q = wv_avx512_exp_f64(_mm512_mul_pd(_mm512_set1_pd(-1.0), _mm512_mul_pd(q, safe_tau)));
    __m512d fwd_S  = _mm512_mul_pd(S, disc_q);
    __m512d fwd_K  = _mm512_mul_pd(K, disc_r);

    __m512d call_pr = _mm512_sub_pd(_mm512_mul_pd(fwd_S, Nd1), _mm512_mul_pd(fwd_K, Nd2));
    __m512d put_pr  = _mm512_sub_pd(_mm512_mul_pd(fwd_K, _mm512_sub_pd(one, Nd2)),
                                    _mm512_mul_pd(fwd_S, _mm512_sub_pd(one, Nd1)));
    __m512d live    = _mm512_mask_blend_pd(call_mask, put_pr, call_pr);
    return _mm512_mask_blend_pd(expired, live, exp_pr);
}

static inline void bsm_implied_vol_vec_x86(const double* restrict S,
                                           const double* restrict K,
                                           const double* restrict tau,
                                           const double* restrict r,
                                           const double* restrict q,
                                           const double* restrict mkt_price,
                                           const uint8_t* restrict is_call,
                                           size_t n,
                                           double* restrict iv_out) {
    WV_REQUIRE(S && K && tau && r && q && mkt_price && is_call && iv_out, "non-null pointers");

    __m512d one = _mm512_set1_pd(1.0);

    #ifdef _OPENMP
    #pragma omp parallel for schedule(static, WV_OMP_CHUNK_SIZE)
    #endif
    for (size_t i = 0; i < (n & ~7ULL); i += 8) {
        __m512d Sv   = _mm512_loadu_pd(S + i);
        __m512d Kv   = _mm512_loadu_pd(K + i);
        __m512d tauv = _mm512_loadu_pd(tau + i);
        __m512d rv   = _mm512_loadu_pd(r + i);
        __m512d qv   = _mm512_loadu_pd(q + i);
        __m512d mv   = _mm512_loadu_pd(mkt_price + i);

        __mmask8 call_mask = 0;
        for (int j = 0; j < 8; j++)
            call_mask |= ((uint8_t)(!!is_call[i+j]) << j);

        double sg[8];
        for (int j = 0; j < 8; j++)
            sg[j] = wv_iv_initial_guess(S[i+j], K[i+j], tau[i+j], r[i+j], q[i+j],
                                        mkt_price[i+j], is_call[i+j]);
        __m512d sigma = _mm512_loadu_pd(sg);

        __m512d safe_tau = _mm512_max_pd(tauv, _mm512_set1_pd(1e-30));
        __m512d sqrt_tau = _mm512_sqrt_pd(safe_tau);
        __m512d logSK    = wv_avx512_log_f64(_mm512_div_pd(Sv, Kv));
        __m512d disc_r   = wv_avx512_exp_f64(_mm512_mul_pd(_mm512_set1_pd(-1.0), _mm512_mul_pd(rv, safe_tau)));
        __m512d disc_q   = wv_avx512_exp_f64(_mm512_mul_pd(_mm512_set1_pd(-1.0), _mm512_mul_pd(qv, safe_tau)));
        __m512d fwd_S    = _mm512_mul_pd(Sv, disc_q);
        __m512d fwd_K    = _mm512_mul_pd(Kv, disc_r);

        __mmask8 active = 0xFF;

        for (int iter = 0; iter < 50; ++iter) {
            __m512d halfv2 = _mm512_mul_pd(_mm512_set1_pd(0.5), _mm512_mul_pd(sigma, sigma));
            __m512d num    = _mm512_fmadd_pd(_mm512_sub_pd(_mm512_add_pd(rv, halfv2), qv),
                                             safe_tau, logSK);
            __m512d d1     = _mm512_div_pd(num, _mm512_mul_pd(sigma, sqrt_tau));
            __m512d d2     = _mm512_sub_pd(d1, _mm512_mul_pd(sigma, sqrt_tau));

            __m512d Nd1   = norm_cdf_avx512(d1);
            __m512d Nd2   = norm_cdf_avx512(d2);
            __m512d call_pr = _mm512_sub_pd(_mm512_mul_pd(fwd_S, Nd1), _mm512_mul_pd(fwd_K, Nd2));
            __m512d put_pr  = _mm512_sub_pd(_mm512_mul_pd(fwd_K, _mm512_sub_pd(one, Nd2)),
                                           _mm512_mul_pd(fwd_S, _mm512_sub_pd(one, Nd1)));
            __m512d price  = _mm512_mask_blend_pd(call_mask, put_pr, call_pr);
            __m512d diff   = _mm512_sub_pd(price, mv);
            __m512d abs_diff = _mm512_abs_pd(diff);

            double ad[8];
            _mm512_storeu_pd(ad, abs_diff);
            int all_cvg = 1;
            for (int j = 0; j < 8; j++) if (ad[j] >= WV_NR_TOL) { all_cvg = 0; break; }
            if (all_cvg) break;

            __m512d vega    = _mm512_mul_pd(fwd_S, _mm512_mul_pd(norm_pdf_avx512(d1), sqrt_tau));
            __mmask8 vega_ok = _mm512_cmp_pd_mask(vega, _mm512_set1_pd(WV_NR_TOL), _CMP_GE_OQ);
            active &= vega_ok;
            if (active == 0) break;

            __m512d step       = _mm512_div_pd(diff, _mm512_max_pd(vega, _mm512_set1_pd(WV_NR_TOL)));
            __m512d half_sigma = _mm512_mul_pd(_mm512_set1_pd(0.5), sigma);
            step = _mm512_min_pd(step, half_sigma);
            step = _mm512_max_pd(step, _mm512_sub_pd(_mm512_setzero_pd(), half_sigma));

            __m512d new_sigma = _mm512_sub_pd(sigma, step);
            new_sigma = _mm512_max_pd(_mm512_set1_pd(0.001), _mm512_min_pd(_mm512_set1_pd(5.0), new_sigma));
            sigma = _mm512_mask_blend_pd(active, sigma, new_sigma);
        }

        _mm512_storeu_pd(iv_out + i, sigma);

        double tv[8];
        _mm512_storeu_pd(tv, tauv);
        for (int j = 0; j < 8; j++)
            if (tv[j] <= 0.0) iv_out[i+j] = 0.0;
    }

    /* Remainder: scalar */
    for (size_t i = (n & ~7ULL); i < n; i++) {
        iv_out[i] = bsm_implied_vol(S[i], K[i], tau[i], r[i], q[i], mkt_price[i], is_call[i]);
    }
}

static inline double implied_correlation_x86(double etf_iv, size_t n_comps,
                                             const ComponentVol* restrict comps) {
    return implied_correlation_scalar(etf_iv, n_comps, comps);
}

static inline void scan_dispersion_opps_x86(
    const double* restrict etf_ivs, const double* restrict etf_vegas, size_t n_etf,
    const ComponentVol* restrict comps, size_t n_comps,
    double target_rho, double thresh, double basket_vega,
    DispersionOpportunity* restrict opps) {
    scan_dispersion_opps_scalar(etf_ivs, etf_vegas, n_etf, comps, n_comps,
                                target_rho, thresh, basket_vega, opps);
}

#endif /* WV_ARCH_* */

/* ===================================================================
 * Arch-neutral public API aliases
 * =================================================================== */
#if defined(WV_ARCH_NEON)
#   define bsm_implied_vol_vec  bsm_implied_vol_vec_arm
#   define scan_dispersion_opps scan_dispersion_opps_arm
#   define implied_correlation  implied_correlation_arm
#elif defined(WV_ARCH_AVX512) || defined(WV_ARCH_AVX2)
#   define bsm_implied_vol_vec  bsm_implied_vol_vec_x86
#   define scan_dispersion_opps scan_dispersion_opps_x86
#   define implied_correlation  implied_correlation_x86
#else
#   define bsm_implied_vol_vec  bsm_implied_vol_vec_scalar
#   define scan_dispersion_opps scan_dispersion_opps_scalar
#   define implied_correlation  implied_correlation_scalar
#endif

#endif /* WONDERVOL_H */
