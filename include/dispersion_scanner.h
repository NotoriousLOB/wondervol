#pragma once
#ifndef DISPERSION_SCANNER_H
#define DISPERSION_SCANNER_H

#include <arm_neon.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <math.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <stdint.h>
#include <errno.h>

#ifndef DISP_STRICT_TRAP
#   define DISP_STRICT_TRAP 0
#endif

#define DISP_REQUIRE(cond, msg) do { \
    if (!(cond)) { \
        static int warned = 0; \
        if (!warned) { \
            fprintf(stderr, "DISPERSION_RIGOUR_VIOLATION: %s (%s:%d) — %s\n", \
                    #cond, __FILE__, __LINE__, msg); \
            warned = 1; \
        } \
        if (DISP_STRICT_TRAP) abort(); \
        errno = EINVAL; \
    } \
} while (0)

/* ===================================================================
 * Constants
 * =================================================================== */
static const double M_INV_SQRT_2PI = 0.39894228040143267794;
static const double AS_P  = 0.2316419;
static const double AS_B1 =  0.319381530;
static const double AS_B2 = -0.356563782;
static const double AS_B3 =  1.781477937;
static const double AS_B4 = -1.821255978;
static const double AS_B5 =  1.330274429;

/* NEON scalar fallback helpers — still shit but correct */
static inline float64x2_t vexpq_f64_custom(float64x2_t x) {
    double v0 = exp(vgetq_lane_f64(x, 0));
    double v1 = exp(vgetq_lane_f64(x, 1));
    return vsetq_lane_f64(v1, vdupq_n_f64(v0), 1);
}

static inline float64x2_t vlogq_f64_custom(float64x2_t x) {
    double v0 = log(vgetq_lane_f64(x, 0));
    double v1 = log(vgetq_lane_f64(x, 1));
    return vsetq_lane_f64(v1, vdupq_n_f64(v0), 1);
}

/* ===================================================================
 * Norm CDF / PDF (NEON)
 * =================================================================== */
static inline float64x2_t norm_cdf_neon(float64x2_t x) {
    float64x2_t zero = vdupq_n_f64(0.0);
    float64x2_t one  = vdupq_n_f64(1.0);
    float64x2_t absx = vabsq_f64(x);

    float64x2_t t = vdivq_f64(one, vfmaq_f64(one, vdupq_n_f64(AS_P), absx));

    float64x2_t t2 = vmulq_f64(t, t);
    float64x2_t t3 = vmulq_f64(t2, t);
    float64x2_t t4 = vmulq_f64(t3, t);
    float64x2_t t5 = vmulq_f64(t4, t);

    float64x2_t poly = vfmaq_f64(vmulq_f64(vdupq_n_f64(AS_B1), t),
                                 vdupq_n_f64(AS_B2), t2);
    poly = vfmaq_f64(poly, vdupq_n_f64(AS_B3), t3);
    poly = vfmaq_f64(poly, vdupq_n_f64(AS_B4), t4);
    poly = vfmaq_f64(poly, vdupq_n_f64(AS_B5), t5);

    float64x2_t pdf = vmulq_f64(vdupq_n_f64(M_INV_SQRT_2PI),
                                vexpq_f64_custom(vmulq_f64(vdupq_n_f64(-0.5), vmulq_f64(absx, absx))));

    float64x2_t tail = vmulq_f64(pdf, poly);
    float64x2_t pos  = vsubq_f64(one, tail);

    uint64x2_t neg_mask = vcltq_f64(x, zero);
    return vbslq_f64(neg_mask, tail, pos);
}

static inline float64x2_t norm_pdf_neon(float64x2_t x) {
    float64x2_t xx = vmulq_f64(x, x);
    return vmulq_f64(vdupq_n_f64(M_INV_SQRT_2PI),
                     vexpq_f64_custom(vmulq_f64(vdupq_n_f64(-0.5), xx)));
}

/* ===================================================================
 * Scalar Norm CDF / PDF (Abramowitz & Stegun 26.2.17)
 * =================================================================== */
static inline double norm_pdf(double x) {
    return M_INV_SQRT_2PI * exp(-0.5 * x * x);
}

static inline double norm_cdf(double x) {
    double absx = fabs(x);
    double t = 1.0 / (1.0 + AS_P * absx);
    double t2 = t * t, t3 = t2 * t, t4 = t3 * t, t5 = t4 * t;
    double poly = AS_B1*t + AS_B2*t2 + AS_B3*t3 + AS_B4*t4 + AS_B5*t5;
    double pdf_val = M_INV_SQRT_2PI * exp(-0.5 * absx * absx);
    double tail = pdf_val * poly;
    return (x < 0.0) ? tail : (1.0 - tail);
}

/* ===================================================================
 * Scalar BSM price (identical to x86 version)
 * =================================================================== */
static inline double bsm_price(double S, double K, double tau, double r, double q,
                               double sigma, int is_call) {
    if (tau <= 0.0 || !isfinite(tau)) {
        return is_call ? fmax(S - K, 0.0) : fmax(K - S, 0.0);
    }
    double d1 = (log(S/K) + (r - q + 0.5*sigma*sigma)*tau) / (sigma * sqrt(tau));
    double d2 = d1 - sigma * sqrt(tau);
    double disc_r = exp(-r * tau);
    double disc_q = exp(-q * tau);
    double Nd1 = norm_cdf(d1);
    double Nd2 = norm_cdf(d2);
    double call = S * disc_q * Nd1 - K * disc_r * Nd2;
    double put  = K * disc_r * (1.0 - Nd2) - S * disc_q * (1.0 - Nd1);
    return is_call ? call : put;
}

/* Scalar IV solver (full) */
static inline double bsm_implied_vol(double S, double K, double tau, double r, double q,
                                     double mkt_price, int is_call) {
    DISP_REQUIRE(S > 0 && K > 0 && isfinite(mkt_price), "invalid BSM params");
    if (tau <= 0.0) return 0.0;

    double intrinsic = is_call ? fmax(S - K, 0.0) : fmax(K - S, 0.0);
    double time_val  = fabs(mkt_price - intrinsic);
    double sigma     = (time_val > 0.0) ? sqrt(2.0 * M_PI / tau) * (time_val / S) : 0.0;
    if (sigma < 0.05) sigma = 0.2; /* B-S guess unreliable for deep OTM */
    sigma = fmax(0.001, fmin(5.0, sigma));

    for (int i = 0; i < 100; ++i) {
        double price = bsm_price(S, K, tau, r, q, sigma, is_call);
        double diff = price - mkt_price;
        if (fabs(diff) < 1e-12) return sigma;
        double d1 = (log(S/K) + (r - q + 0.5*sigma*sigma)*tau) / (sigma * sqrt(tau));
        double vega = S * exp(-q*tau) * norm_pdf(d1) * sqrt(tau);
        if (vega < 1e-12) break;
        double step = diff / vega;
        /* Dampen: limit step to half of current sigma to prevent oscillation */
        if (step > 0.5 * sigma) step = 0.5 * sigma;
        if (step < -0.5 * sigma) step = -0.5 * sigma;
        sigma -= step;
        sigma = fmax(0.001, fmin(5.0, sigma));
        if (fabs(diff) < 1e-9) break;
    }
    return sigma;
}

/* ===================================================================
 * Vectorized BSM price for NEON (2-wide)
 * =================================================================== */
static inline float64x2_t bsm_price_neon(float64x2_t S, float64x2_t K, float64x2_t tau,
                                         float64x2_t r, float64x2_t q, float64x2_t sigma,
                                         uint64x2_t call_mask) {
    float64x2_t zero = vdupq_n_f64(0.0);
    float64x2_t one  = vdupq_n_f64(1.0);

    /* Expired intrinsic values */
    float64x2_t call_int = vmaxq_f64(vsubq_f64(S, K), zero);
    float64x2_t put_int  = vmaxq_f64(vsubq_f64(K, S), zero);
    float64x2_t expired_price = vbslq_f64(call_mask, call_int, put_int);

    uint64x2_t expired = vcleq_f64(tau, zero);

    /* Fast path: both lanes expired */
    if (vgetq_lane_u64(expired, 0) && vgetq_lane_u64(expired, 1)) {
        return expired_price;
    }

    /* BSM math for live options — guard tau to avoid sqrt(0) */
    float64x2_t safe_tau = vmaxq_f64(tau, vdupq_n_f64(1e-30));
    float64x2_t sqrt_tau = vsqrtq_f64(safe_tau);
    float64x2_t sig_sqrt_tau = vmulq_f64(sigma, sqrt_tau);
    float64x2_t logSK = vlogq_f64_custom(vdivq_f64(S, K));
    float64x2_t half_sig2 = vmulq_f64(vdupq_n_f64(0.5), vmulq_f64(sigma, sigma));
    float64x2_t drift = vfmaq_f64(logSK, vsubq_f64(vaddq_f64(r, half_sig2), q), safe_tau);
    float64x2_t d1 = vdivq_f64(drift, sig_sqrt_tau);
    float64x2_t d2 = vsubq_f64(d1, sig_sqrt_tau);

    float64x2_t Nd1 = norm_cdf_neon(d1);
    float64x2_t Nd2 = norm_cdf_neon(d2);

    float64x2_t disc_r = vexpq_f64_custom(vnegq_f64(vmulq_f64(r, safe_tau)));
    float64x2_t disc_q = vexpq_f64_custom(vnegq_f64(vmulq_f64(q, safe_tau)));

    float64x2_t fwd_S = vmulq_f64(S, disc_q);
    float64x2_t fwd_K = vmulq_f64(K, disc_r);

    float64x2_t call_pr = vsubq_f64(vmulq_f64(fwd_S, Nd1), vmulq_f64(fwd_K, Nd2));
    float64x2_t put_pr  = vsubq_f64(vmulq_f64(fwd_K, vsubq_f64(one, Nd2)),
                                     vmulq_f64(fwd_S, vsubq_f64(one, Nd1)));
    float64x2_t live_price = vbslq_f64(call_mask, call_pr, put_pr);

    /* Blend: expired lanes get intrinsic, live lanes get BSM price */
    return vbslq_f64(expired, expired_price, live_price);
}

/* ===================================================================
 * NEON Implied Vol Solver — production ready
 * =================================================================== */
#ifndef DISP_OMP_CHUNK_SIZE
#define DISP_OMP_CHUNK_SIZE 256
#endif

static inline void bsm_implied_vol_vec_arm(const double* restrict S,
                                           const double* restrict K,
                                           const double* restrict tau,
                                           const double* restrict r,
                                           const double* restrict q,
                                           const double* restrict mkt_price,
                                           const uint8_t* restrict is_call,
                                           size_t n,
                                           double* restrict iv_out) {
    DISP_REQUIRE(S && K && tau && r && q && mkt_price && is_call && iv_out, "non-null pointers");

    #ifdef _OPENMP
    #pragma omp parallel for schedule(static, DISP_OMP_CHUNK_SIZE)
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

        float64x2_t zero = vdupq_n_f64(0.0);
        float64x2_t one  = vdupq_n_f64(1.0);
        float64x2_t intrinsic = vbslq_f64(call_mask,
                                          vmaxq_f64(vsubq_f64(Sv, Kv), zero),
                                          vmaxq_f64(vsubq_f64(Kv, Sv), zero));

        /* A5: Use vabsq for correct abs; guard tau against zero/negative */
        float64x2_t time_val = vabsq_f64(vsubq_f64(mv, intrinsic));
        float64x2_t safe_tau = vmaxq_f64(tauv, vdupq_n_f64(1e-30));
        float64x2_t raw_sigma = vsqrtq_f64(vdivq_f64(vmulq_f64(vdupq_n_f64(2.0 * M_PI), time_val),
                                                      vmulq_f64(Sv, safe_tau)));
        /* NaN/low fallback: substitute 0.2 for non-finite or unreliably small guesses */
        uint64x2_t is_finite = vandq_u64(vcleq_f64(raw_sigma, vdupq_n_f64(1e30)),
                                         vcgeq_f64(raw_sigma, zero));
        float64x2_t sigma = vbslq_f64(is_finite, raw_sigma, vdupq_n_f64(0.2));
        uint64x2_t too_low = vcltq_f64(sigma, vdupq_n_f64(0.05));
        sigma = vbslq_f64(too_low, vdupq_n_f64(0.2), sigma);
        sigma = vmaxq_f64(vdupq_n_f64(0.001), vminq_f64(vdupq_n_f64(5.0), sigma));

        float64x2_t sqrt_tau = vsqrtq_f64(safe_tau);

        /* B1: Hoist loop-invariant transcendentals */
        float64x2_t logSK = vlogq_f64_custom(vdivq_f64(Sv, Kv));
        float64x2_t disc_r = vexpq_f64_custom(vnegq_f64(vmulq_f64(rv, safe_tau)));
        float64x2_t disc_q = vexpq_f64_custom(vnegq_f64(vmulq_f64(qv, safe_tau)));
        float64x2_t fwd_S = vmulq_f64(Sv, disc_q);
        float64x2_t fwd_K = vmulq_f64(Kv, disc_r);

        for (int iter = 0; iter < 50; ++iter) {
            float64x2_t halfv2 = vmulq_f64(vdupq_n_f64(0.5), vmulq_f64(sigma, sigma));
            float64x2_t num = vfmaq_f64(logSK, vsubq_f64(vaddq_f64(rv, halfv2), qv), safe_tau);
            float64x2_t d1 = vdivq_f64(num, vmulq_f64(sigma, sqrt_tau));
            float64x2_t d2 = vsubq_f64(d1, vmulq_f64(sigma, sqrt_tau));

            float64x2_t Nd1 = norm_cdf_neon(d1);
            float64x2_t Nd2 = norm_cdf_neon(d2);

            float64x2_t call_pr = vsubq_f64(vmulq_f64(fwd_S, Nd1), vmulq_f64(fwd_K, Nd2));
            float64x2_t put_pr  = vsubq_f64(vmulq_f64(fwd_K, vsubq_f64(one, Nd2)),
                                            vmulq_f64(fwd_S, vsubq_f64(one, Nd1)));
            float64x2_t price = vbslq_f64(call_mask, call_pr, put_pr);

            float64x2_t diff = vsubq_f64(price, mv);
            float64x2_t abs_diff = vabsq_f64(diff);

            /* early convergence exit */
            if (vgetq_lane_f64(abs_diff, 0) < 1e-12 && vgetq_lane_f64(abs_diff, 1) < 1e-12) {
                break;
            }

            float64x2_t vega = vmulq_f64(fwd_S, vmulq_f64(norm_pdf_neon(d1), sqrt_tau));
            float64x2_t vega_clamp = vmaxq_f64(vega, vdupq_n_f64(1e-12));

            float64x2_t step = vdivq_f64(diff, vega_clamp);
            /* Dampen: limit step to half of current sigma to prevent oscillation */
            float64x2_t half_sigma = vmulq_f64(vdupq_n_f64(0.5), sigma);
            step = vminq_f64(step, half_sigma);
            step = vmaxq_f64(step, vnegq_f64(half_sigma));
            sigma = vsubq_f64(sigma, step);
            sigma = vmaxq_f64(vdupq_n_f64(0.001), vminq_f64(vdupq_n_f64(5.0), sigma));

            /* convergence check */
            if (vgetq_lane_f64(abs_diff, 0) < 1e-9 && vgetq_lane_f64(abs_diff, 1) < 1e-9) break;
        }

        vst1q_f64(iv_out + i, sigma);

        /* A6: Fix up expired lanes — IV = 0.0 when tau <= 0 */
        if (vgetq_lane_f64(tauv, 0) <= 0.0) iv_out[i] = 0.0;
        if (vgetq_lane_f64(tauv, 1) <= 0.0) iv_out[i + 1] = 0.0;
    }

    /* scalar tail */
    if (n % 2 != 0) {
        size_t last = n - 1;
        iv_out[last] = bsm_implied_vol(S[last], K[last], tau[last], r[last], q[last],
                                       mkt_price[last], is_call[last]);
    }
}

/* ===================================================================
 * Dispersion scanner (identical math to x86 version)
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

static inline double implied_correlation_arm(double etf_iv, size_t n_comps,
                                             const ComponentVol* restrict comps) {
    DISP_REQUIRE(etf_iv > 0.0 && comps, "valid args");
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

static inline void scan_dispersion_opps_arm(
    const double* restrict etf_ivs, const double* restrict etf_vegas, size_t n_etf,
    const ComponentVol* restrict comps, size_t n_comps,
    double target_rho, double thresh, double basket_vega,
    DispersionOpportunity* restrict opps) {

    DISP_REQUIRE(etf_ivs && etf_vegas && comps && opps && n_etf > 0, "valid args");

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

#endif /* DISPERSION_SCANNER_H */
