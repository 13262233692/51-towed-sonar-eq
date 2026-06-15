#define _USE_MATH_DEFINES
#include "doppler_resampler.h"
#include <cmath>
#include <algorithm>
#include <numeric>
#include <chrono>
#include <stdexcept>

namespace sonar {

DopplerResampler::DopplerResampler()
    : sample_rate_hz_(48000.0), farrow_order_(3) {}

void DopplerResampler::set_sample_rate(double sample_rate_hz) {
    sample_rate_hz_ = sample_rate_hz;
}

void DopplerResampler::set_farrow_order(size_t order) {
    farrow_order_ = std::min((size_t)5, std::max((size_t)1, order));
}

double DopplerResampler::lagrange_kernel(double mu, size_t order, size_t tap) {
    double result = 1.0;
    int center = (int)(FIR_HALF_TAPS);
    for (size_t k = 0; k <= order; ++k) {
        if (k == tap) continue;
        double n_k = (double)((int)k - center);
        double n_t = (double)((int)tap - center);
        result *= (mu - n_k) / (n_t - n_k);
    }
    return result;
}

std::vector<std::vector<double>> DopplerResampler::compute_farrow_coeff_table(
    double rate_ratio) {
    size_t M = FIR_HALF_TAPS * 2 + 1;
    size_t N = farrow_order_ + 1;
    std::vector<std::vector<double>> table(N, std::vector<double>(M, 0.0));
    std::vector<double> mus(N);
    for (size_t i = 0; i < N; ++i) mus[i] = (double)i / (double)(N - 1);
    for (size_t n = 0; n < N; ++n) {
        double mu = mus[n] * rate_ratio;
        for (size_t tap = 0; tap < M; ++tap) {
            double sum = 0.0;
            for (size_t k = 0; k <= farrow_order_; ++k) {
                double coeff_k = lagrange_kernel(mus[k], farrow_order_, tap);
                double lag_basis = 1.0;
                for (size_t j = 0; j <= farrow_order_; ++j) {
                    if (j == k) continue;
                    lag_basis *= (mu - mus[j]) / (mus[k] - mus[j]);
                }
                sum += coeff_k * lag_basis;
            }
            table[n][tap] = sum;
        }
    }
    return table;
}

std::complex<double> DopplerResampler::farrow_interpolate(
    const std::vector<std::complex<double>>& buf,
    double mu,
    const std::vector<std::vector<double>>& coeffs) {
    std::complex<double> result(0, 0);
    size_t M = FIR_HALF_TAPS * 2 + 1;
    for (size_t tap = 0; tap < M; ++tap) {
        std::complex<double> tap_val(0, 0);
        for (int n = (int)farrow_order_; n >= 0; --n) {
            size_t idx = (size_t)n;
            tap_val = tap_val * mu + std::complex<double>(coeffs[idx][tap], 0.0);
        }
        SONAR_BOUNDS_CHECK(tap, buf.size());
        result += tap_val * buf[tap];
    }
    return result;
}

std::vector<ResamplingPoly> DopplerResampler::design_resampling_polynomials(
    const std::vector<double>& instantaneous_doppler_hz,
    double carrier_frequency_hz,
    size_t signal_length,
    double poly_update_interval) {
    std::vector<ResamplingPoly> polys;
    if (instantaneous_doppler_hz.empty() || signal_length == 0) return polys;
    size_t num_polys = (size_t)std::ceil((double)signal_length / poly_update_interval);
    polys.reserve(num_polys);
    for (size_t p = 0; p < num_polys; ++p) {
        ResamplingPoly poly;
        poly.valid_start_sample = (double)p * poly_update_interval;
        poly.valid_end_sample = std::min((double)(p + 1) * poly_update_interval, (double)signal_length);
        size_t idx_start = (size_t)std::min((size_t)p * (size_t)poly_update_interval,
                                            instantaneous_doppler_hz.size() - 1);
        size_t idx_end = (size_t)std::min(idx_start + (size_t)poly_update_interval,
                                          instantaneous_doppler_hz.size() - 1);
        double avg_doppler = 0.0;
        size_t cnt = 0;
        for (size_t i = idx_start; i <= idx_end && i < instantaneous_doppler_hz.size(); ++i) {
            avg_doppler += instantaneous_doppler_hz[i];
            cnt++;
        }
        if (cnt > 0) avg_doppler /= (double)cnt;
        double fc = std::max(1.0, carrier_frequency_hz);
        double rate_ratio = fc / (fc + avg_doppler);
        poly.rate_ratio = rate_ratio;
        double t0 = poly.valid_start_sample;
        double t1 = poly.valid_end_sample;
        double d0 = (idx_start < instantaneous_doppler_hz.size()) ?
                    instantaneous_doppler_hz[idx_start] : 0.0;
        double d1 = (idx_end < instantaneous_doppler_hz.size()) ?
                    instantaneous_doppler_hz[idx_end] : 0.0;
        double s0 = fc / std::max(1e-6, fc + d0);
        double s1 = fc / std::max(1e-6, fc + d1);
        double dt = std::max(1e-6, t1 - t0);
        double c0 = s0;
        double c1 = (s1 - s0) / dt;
        double c2 = 0.0;
        double c3 = 0.0;
        if (farrow_order_ >= 3) {
            size_t idx_mid = (idx_start + idx_end) / 2;
            double dm = (idx_mid < instantaneous_doppler_hz.size()) ?
                        instantaneous_doppler_hz[idx_mid] : 0.0;
            double sm = fc / std::max(1e-6, fc + dm);
            double tm = (t0 + t1) / 2.0;
            double denom = (tm - t0) * (tm - t1) * (t0 - t1);
            if (std::abs(denom) > 1e-12) {
                c2 = (sm * (t1 - t0) - c0 * (tm - t1) - c1 * (tm - t0) * (t1 - t0)) / denom;
                c3 = 0.0;
            }
        }
        poly.coeffs[0] = c0;
        poly.coeffs[1] = c1;
        if (FARROW_POLY_ORDER >= 2) poly.coeffs[2] = c2;
        if (FARROW_POLY_ORDER >= 3) poly.coeffs[3] = c3;
        for (size_t k = FARROW_POLY_ORDER + 1; k < poly.coeffs.size(); ++k) poly.coeffs[k] = 0.0;
        poly.fractional_delay = 0.0;
        polys.push_back(poly);
    }
    return polys;
}

std::vector<double> DopplerResampler::estimate_instantaneous_doppler(
    const std::vector<AcousticRay>& rays,
    size_t signal_length,
    double sample_rate_hz) {
    std::vector<double> dopplers(signal_length, 0.0);
    if (rays.empty()) return dopplers;
    double total_weight = 0.0;
    double weighted_doppler = 0.0;
    for (const auto& r : rays) {
        double w = std::abs(r.complex_weight);
        if (!std::isfinite(w)) continue;
        total_weight += w;
        weighted_doppler += w * r.doppler_shift_hz;
    }
    double mean_dop = (total_weight > 1e-12) ? weighted_doppler / total_weight : 0.0;
    double max_delta = 0.0;
    for (const auto& r : rays) {
        if (std::isfinite(r.doppler_shift_hz)) {
            max_delta = std::max(max_delta, std::abs(r.doppler_shift_hz - mean_dop));
        }
    }
    for (size_t n = 0; n < signal_length; ++n) {
        double t = (double)n / std::max(1.0, sample_rate_hz);
        double mod = std::sin(2.0 * M_PI * 2.5 * t) * 0.15 + 1.0;
        double extra = 0.0;
        if (rays.size() > 1) {
            size_t idx = n % rays.size();
            extra = (rays[idx].doppler_shift_hz - mean_dop) * 0.2;
        }
        dopplers[n] = mean_dop * mod + extra;
    }
    return dopplers;
}

std::vector<std::complex<double>> DopplerResampler::apply_farrow_resampling(
    const std::vector<std::complex<double>>& input,
    const std::vector<ResamplingPoly>& polys) {
    if (input.empty() || polys.empty()) return input;
    size_t N_in = input.size();
    double avg_rate = 0.0;
    for (const auto& p : polys) avg_rate += p.rate_ratio;
    avg_rate /= (double)polys.size();
    size_t N_out = (size_t)std::ceil((double)N_in / std::max(1e-6, avg_rate)) + 10;
    std::vector<std::complex<double>> output;
    output.reserve(N_out);
    double t_out = 0.0;
    size_t half_tap = FIR_HALF_TAPS;
    auto eval_poly = [&](const ResamplingPoly& p, double t_sample) -> double {
        double dt = t_sample - p.valid_start_sample;
        double result = 0.0;
        for (int k = (int)FARROW_POLY_ORDER; k >= 0; --k) {
            result = result * dt + p.coeffs[k];
        }
        return result;
    };
    auto find_poly = [&](double t) -> const ResamplingPoly& {
        for (const auto& p : polys) {
            if (t >= p.valid_start_sample && t < p.valid_end_sample) return p;
        }
        return polys.back();
    };
    size_t M_total = FIR_HALF_TAPS * 2 + 1;
    size_t max_iters = N_in * 3 + 1000;
    for (size_t iter = 0; iter < max_iters; ++iter) {
        if (t_out >= (double)N_in - (double)half_tap - 1) break;
        double floor_t = std::floor(t_out);
        size_t base = (size_t)floor_t;
        double mu = t_out - floor_t;
        const ResamplingPoly& poly = find_poly(t_out);
        double cur_rate = eval_poly(poly, t_out);
        cur_rate = std::max(0.1, std::min(10.0, cur_rate));
        std::vector<std::complex<double>> buf(M_total, std::complex<double>(0, 0));
        for (size_t k = 0; k < M_total; ++k) {
            long long idx = (long long)base + (long long)k - (long long)half_tap;
            if (idx < 0 || idx >= (long long)N_in) {
                buf[k] = std::complex<double>(0.0, 0.0);
            } else {
                buf[k] = input[(size_t)idx];
            }
        }
        auto table = compute_farrow_coeff_table(cur_rate);
        output.push_back(farrow_interpolate(buf, mu, table));
        t_out += cur_rate;
    }
    return output;
}

ComplexVector DopplerResampler::phase_derotate(
    const ComplexVector& signal,
    const std::vector<double>& doppler_hz,
    double sample_rate_hz) {
    ComplexVector out(signal.size());
    double phase = 0.0;
    for (size_t n = 0; n < signal.size(); ++n) {
        double d = (n < doppler_hz.size()) ? doppler_hz[n] : 0.0;
        double dphi = 2.0 * M_PI * d / std::max(1.0, sample_rate_hz);
        phase += dphi;
        phase = std::fmod(phase, 2.0 * M_PI);
        std::complex<double> rot(std::cos(-phase), std::sin(-phase));
        out[n] = signal[n] * rot;
    }
    return out;
}

DopplerCompensationResult DopplerResampler::compensate(
    const std::vector<std::complex<double>>& input_signal,
    const PlatformState& state,
    const ArrayDeformation& deformation,
    AcousticRayTracer& tracer,
    size_t num_rays) {
    auto t0 = std::chrono::high_resolution_clock::now();
    DopplerCompensationResult result;
    if (input_signal.empty()) {
        result.processing_time_ms = 0.0;
        result.total_compensation_gain_db = 0.0;
        result.residual_drift_ppm = 0.0;
        return result;
    }
    result.active_rays = tracer.trace_rays(state, deformation, num_rays / 4, 2);
    result.instantaneous_doppler_hz = estimate_instantaneous_doppler(
        result.active_rays, input_signal.size(), state.baseband_sample_rate_hz);
    result.radial_velocity_track_mps.resize(input_signal.size());
    for (size_t n = 0; n < input_signal.size(); ++n) {
        result.radial_velocity_track_mps[n] =
            result.instantaneous_doppler_hz[n] / std::max(1.0, state.carrier_frequency_hz)
            * SPEED_OF_SOUND_WATER;
    }
    result.poly_series = design_resampling_polynomials(
        result.instantaneous_doppler_hz, state.carrier_frequency_hz, input_signal.size(), 128.0);
    ComplexVector derotated = phase_derotate(input_signal,
        result.instantaneous_doppler_hz, state.baseband_sample_rate_hz);
    result.resampled_signal = apply_farrow_resampling(derotated, result.poly_series);
    double mean_dop = 0.0;
    for (double d : result.instantaneous_doppler_hz) mean_dop += d;
    mean_dop /= (double)result.instantaneous_doppler_hz.size();
    result.total_compensation_gain_db = 20.0 * std::log10(
        std::max(1e-6, 1.0 + std::abs(mean_dop) / std::max(1.0, state.carrier_frequency_hz)));
    size_t N = result.instantaneous_doppler_hz.size();
    double residual = 0.0;
    for (size_t i = N - std::min((size_t)100, N); i < N; ++i) {
        residual += result.instantaneous_doppler_hz[i] * result.instantaneous_doppler_hz[i];
    }
    residual = std::sqrt(residual / 100.0);
    result.residual_drift_ppm = residual / std::max(1.0, state.carrier_frequency_hz) * 1e6;
    auto t1 = std::chrono::high_resolution_clock::now();
    result.processing_time_ms =
        std::chrono::duration<double, std::milli>(t1 - t0).count();
    return result;
}

}
