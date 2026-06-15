#include "cma_equalizer.h"
#include "sonar_types.h"
#include <chrono>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <numeric>
#include <cstring>
#include <limits>

namespace sonar {

namespace {

BFPBlock compute_bfp_block(const ComplexVector& signal,
                           size_t block_start, size_t block_len) {
    BFPBlock blk;
    blk.mantissa.resize(block_len, std::complex<double>(0.0, 0.0));

    double max_abs = 0.0;
    for (size_t i = 0; i < block_len; ++i) {
        const size_t abs_idx = block_start + i;
        SONAR_BOUNDS_CHECK(abs_idx, signal.size());
        double a = std::abs(signal[abs_idx]);
        if (a > max_abs) max_abs = a;
    }

    if (max_abs < std::numeric_limits<double>::min()) {
        blk.exponent = 0;
        return blk;
    }

    blk.exponent = 0;
    double s = max_abs;
    while (s >= 2.0) {
        s *= 0.5;
        blk.exponent++;
    }
    while (s < 0.5 && blk.exponent > -1022) {
        s *= 2.0;
        blk.exponent--;
    }

    const double inv_scale = std::ldexp(1.0, -blk.exponent);
    for (size_t i = 0; i < block_len; ++i) {
        const size_t sig_idx = block_start + i;
        SONAR_BOUNDS_CHECK(sig_idx, signal.size());
        SONAR_BOUNDS_CHECK(i, blk.mantissa.size());
        blk.mantissa[i] = signal[sig_idx] * inv_scale;
    }

    return blk;
}

ComplexVector reconstruct_bfp_signal(const std::vector<BFPBlock>& blocks,
                                     size_t total_samples) {
    ComplexVector out(total_samples, std::complex<double>(0.0, 0.0));
    size_t cursor = 0;

    for (const auto& blk : blocks) {
        const double scale = std::ldexp(1.0, blk.exponent);
        const size_t blk_len = blk.mantissa.size();

        SONAR_ALIGNED_BLOCK_LOOP(k, blk_len, BFP_BLOCK_SIZE) {
            for (size_t j = 0; j < BFP_BLOCK_SIZE; ++j) {
                const size_t gi = cursor + k + j;
                SONAR_BOUNDS_CHECK(gi, out.size());
                SONAR_BOUNDS_CHECK(k + j, blk.mantissa.size());
                out[gi] = blk.mantissa[k + j] * scale;
            }
        }

        SONAR_ALIGNED_REMAINDER_LOOP(k, blk_len, BFP_BLOCK_SIZE) {
            const size_t gi = cursor + k;
            SONAR_BOUNDS_CHECK(gi, out.size());
            SONAR_BOUNDS_CHECK(k, blk.mantissa.size());
            out[gi] = blk.mantissa[k] * scale;
        }

        cursor += blk_len;
    }

    return out;
}

void shift_signal_buffer(ComplexVector& buf, size_t M) {
    if (M < 2) return;
    SONAR_BOUNDS_CHECK(M - 1, buf.size());

    const size_t aligned_end = ((M - 1) / BFP_BLOCK_SIZE) * BFP_BLOCK_SIZE;

    for (size_t k = M - 1; k > aligned_end && k > 0; --k) {
        SONAR_BOUNDS_CHECK(k, buf.size());
        SONAR_BOUNDS_CHECK(k - 1, buf.size());
        buf[k] = buf[k - 1];
    }

    for (size_t block_base = aligned_end;
         block_base + BFP_BLOCK_SIZE <= M - 1 && block_base >= BFP_BLOCK_SIZE;
         block_base -= BFP_BLOCK_SIZE) {
        for (size_t j = 0; j < BFP_BLOCK_SIZE; ++j) {
            const size_t dst = block_base + BFP_BLOCK_SIZE - 1 - j;
            const size_t src = dst - 1;
            SONAR_BOUNDS_CHECK(dst, buf.size());
            SONAR_BOUNDS_CHECK(src, buf.size());
            buf[dst] = buf[src];
        }
    }
}

}

CMAEqualizer::CMAEqualizer(const CMAConfig& config)
    : config_(config) {
    ray_tracer_ = std::make_unique<AcousticRayTracer>();
    resampler_ = std::make_unique<DopplerResampler>();
}

void CMAEqualizer::set_config(const CMAConfig& config) {
    config_ = config;
}

const CMAConfig& CMAEqualizer::get_config() const {
    return config_;
}

void CMAEqualizer::set_platform_state(const PlatformState& state) {
    cached_platform_ = state;
    resampler_->set_sample_rate(state.baseband_sample_rate_hz);
}

void CMAEqualizer::set_array_deformation(const ArrayDeformation& deform) {
    cached_deformation_ = deform;
}

DopplerCompensationResult CMAEqualizer::get_last_doppler_compensation() const {
    return last_doppler_result_;
}

std::complex<double> CMAEqualizer::fir_convolve(
    const ComplexVector& weights,
    const ComplexVector& signal,
    size_t index) const {

    std::complex<double> sum(0.0, 0.0);
    const size_t M = weights.size();

    SONAR_BOUNDS_CHECK(index, signal.size());

    SONAR_ALIGNED_BLOCK_LOOP(k, std::min(M, index + 1), BFP_BLOCK_SIZE) {
        for (size_t j = 0; j < BFP_BLOCK_SIZE; ++j) {
            const size_t kk = k + j;
            if (kk >= M || kk > index) break;
            SONAR_BOUNDS_CHECK(kk, weights.size());
            SONAR_BOUNDS_CHECK(index - kk, signal.size());
            sum += weights[kk] * signal[index - kk];
        }
    }

    SONAR_ALIGNED_REMAINDER_LOOP(k, std::min(M, index + 1), BFP_BLOCK_SIZE) {
        if (k >= M || k > index) continue;
        SONAR_BOUNDS_CHECK(k, weights.size());
        SONAR_BOUNDS_CHECK(index - k, signal.size());
        sum += weights[k] * signal[index - k];
    }

    return sum;
}

void CMAEqualizer::update_weights(
    ComplexVector& weights,
    std::complex<double> error,
    const ComplexVector& signal_buffer,
    double step_size,
    double doppler_bias) const {

    const size_t M = weights.size();
    const double leakage = config_.leakage;

    if (signal_buffer.size() < M) {
        throw std::out_of_range("CMAEqualizer::update_weights: signal_buffer too small");
    }

    double phase_rot = 0.0;
    if (std::abs(doppler_bias) > 1e-12) {
        phase_rot = doppler_bias * config_.doppler_feedback_gain;
    }

    auto apply_weight_update = [&](size_t kk) {
        SONAR_BOUNDS_CHECK(kk, weights.size());
        SONAR_BOUNDS_CHECK(kk, signal_buffer.size());
        std::complex<double> grad = error * std::conj(signal_buffer[kk]);

        const double grad_clip = 1.0 / std::max(step_size, 1e-10);
        double gnorm = std::abs(grad);
        if (gnorm > grad_clip && gnorm > 0) {
            grad *= (grad_clip / gnorm);
        }

        std::complex<double> w_new =
            (1.0 - leakage) * weights[kk] - step_size * grad;

        if (std::abs(phase_rot) > 1e-12 && kk > 0 && kk + 1 < M) {
            size_t prev_idx = kk - 1;
            size_t next_idx = kk + 1;
            SONAR_BOUNDS_CHECK(prev_idx, weights.size());
            SONAR_BOUNDS_CHECK(next_idx, weights.size());
            double interpolate = phase_rot * 0.5;
            auto blended = weights[prev_idx] * (0.5 - interpolate)
                         + weights[kk] * (1.0 - std::abs(interpolate) * 0.3)
                         + weights[next_idx] * (0.5 + interpolate);
            w_new = w_new * 0.7 + blended * 0.3;
        }

        if (!std::isfinite(w_new.real()) || !std::isfinite(w_new.imag())) {
            w_new = std::complex<double>(0.0, 0.0);
        }

        double wmag = std::abs(w_new);
        const double w_clip = 100.0;
        if (wmag > w_clip && wmag > 0) {
            w_new *= (w_clip / wmag);
        }

        weights[kk] = w_new;
    };

    SONAR_ALIGNED_BLOCK_LOOP(k, M, BFP_BLOCK_SIZE) {
        for (size_t j = 0; j < BFP_BLOCK_SIZE; ++j) {
            apply_weight_update(k + j);
        }
    }

    SONAR_ALIGNED_REMAINDER_LOOP(k, M, BFP_BLOCK_SIZE) {
        apply_weight_update(k);
    }
}

ComplexVector CMAEqualizer::apply_filter(
    const ComplexVector& signal,
    const ComplexVector& weights) const {

    const size_t N = signal.size();
    ComplexVector output(N, std::complex<double>(0.0, 0.0));

    for (size_t n = 0; n < N; ++n) {
        output[n] = fir_convolve(weights, signal, n);
    }

    return output;
}

CMAEqualizerResult CMAEqualizer::equalize(const ComplexVector& received_signal) {
    auto t0 = std::chrono::high_resolution_clock::now();

    if (received_signal.empty()) {
        throw std::invalid_argument("Received signal is empty");
    }

    const size_t N = received_signal.size();
    const size_t M = config_.filter_taps;
    if (M == 0) {
        throw std::invalid_argument("filter_taps must be > 0");
    }

    const double R2 = config_.modulus * config_.modulus;
    const double mu = config_.step_size;
    const bool use_bfp = config_.enable_bfp_normalization;

    ComplexVector working_signal;
    std::vector<BFPBlock> bfp_blocks;

    if (use_bfp) {
        bfp_blocks.reserve((N + BFP_BLOCK_SIZE - 1) / BFP_BLOCK_SIZE);
        for (size_t bs = 0; bs < N; bs += BFP_BLOCK_SIZE) {
            const size_t blen = std::min(BFP_BLOCK_SIZE, N - bs);
            bfp_blocks.push_back(compute_bfp_block(received_signal, bs, blen));
        }
        working_signal = reconstruct_bfp_signal(bfp_blocks, N);
    } else {
        working_signal = received_signal;
    }

    ComplexVector weights(M, std::complex<double>(0.0, 0.0));
    const size_t center_tap = M / 2;
    SONAR_BOUNDS_CHECK(center_tap, weights.size());
    weights[center_tap] = std::complex<double>(1.0, 0.0);

    ComplexVector equalized(N, std::complex<double>(0.0, 0.0));
    std::vector<double> convergence_curve;
    convergence_curve.reserve(config_.max_iterations);

    bool converged = false;
    size_t iter = 0;
    double prev_mse = 1e300;
    double final_mse = 0.0;

    ComplexVector signal_buffer(M, std::complex<double>(0.0, 0.0));

    for (; iter < config_.max_iterations; ++iter) {
        double epoch_mse = 0.0;
        size_t count = 0;

        if (use_bfp && (iter % 10 == 0) && iter > 0) {
            for (size_t bs = 0; bs < N; bs += BFP_BLOCK_SIZE) {
                const size_t blen = std::min(BFP_BLOCK_SIZE, N - bs);
                bfp_blocks[bs / BFP_BLOCK_SIZE] =
                    compute_bfp_block(working_signal, bs, blen);
            }
            working_signal = reconstruct_bfp_signal(bfp_blocks, N);
        }

        std::fill(signal_buffer.begin(), signal_buffer.end(),
                  std::complex<double>(0.0, 0.0));

        for (size_t n = 0; n < N; ++n) {
            for (size_t k = M - 1; k > 0; --k) {
                SONAR_BOUNDS_CHECK(k, signal_buffer.size());
                SONAR_BOUNDS_CHECK(k - 1, signal_buffer.size());
                signal_buffer[k] = signal_buffer[k - 1];
            }
            SONAR_BOUNDS_CHECK(0, signal_buffer.size());
            SONAR_BOUNDS_CHECK(n, working_signal.size());
            signal_buffer[0] = working_signal[n];

            std::complex<double> y(0.0, 0.0);
            SONAR_ALIGNED_BLOCK_LOOP(k, M, BFP_BLOCK_SIZE) {
                for (size_t j = 0; j < BFP_BLOCK_SIZE; ++j) {
                    const size_t kk = k + j;
                    SONAR_BOUNDS_CHECK(kk, weights.size());
                    SONAR_BOUNDS_CHECK(kk, signal_buffer.size());
                    y += weights[kk] * signal_buffer[kk];
                }
            }
            SONAR_ALIGNED_REMAINDER_LOOP(k, M, BFP_BLOCK_SIZE) {
                SONAR_BOUNDS_CHECK(k, weights.size());
                SONAR_BOUNDS_CHECK(k, signal_buffer.size());
                y += weights[k] * signal_buffer[k];
            }

            SONAR_BOUNDS_CHECK(n, equalized.size());
            equalized[n] = y;

            double mag_sq = std::norm(y);
            if (!std::isfinite(mag_sq) || mag_sq > 1e12) {
                for (size_t kk = 0; kk < M; ++kk) {
                    SONAR_BOUNDS_CHECK(kk, weights.size());
                    weights[kk] = std::complex<double>(0.0, 0.0);
                }
                const size_t center = M / 2;
                SONAR_BOUNDS_CHECK(center, weights.size());
                weights[center] = std::complex<double>(1.0, 0.0);
                std::fill(signal_buffer.begin(), signal_buffer.end(),
                          std::complex<double>(0.0, 0.0));
                continue;
            }
            std::complex<double> e = y * (mag_sq - R2);

            epoch_mse += std::norm(e);
            ++count;

            update_weights(weights, e, signal_buffer, mu, 0.0);
        }

        if (count > 0) {
            epoch_mse /= count;
        }
        convergence_curve.push_back(epoch_mse);
        final_mse = epoch_mse;

        double delta = std::abs(prev_mse - epoch_mse) / (prev_mse + 1e-30);
        if (iter > 10 && delta < config_.convergence_threshold) {
            converged = true;
            break;
        }
        prev_mse = epoch_mse;
    }

    for (size_t n = 0; n < N; ++n) {
        std::complex<double> y(0.0, 0.0);
        SONAR_ALIGNED_BLOCK_LOOP(k, std::min(M, n + 1), BFP_BLOCK_SIZE) {
            for (size_t j = 0; j < BFP_BLOCK_SIZE; ++j) {
                const size_t kk = k + j;
                if (kk >= M || kk > n) break;
                SONAR_BOUNDS_CHECK(kk, weights.size());
                SONAR_BOUNDS_CHECK(n - kk, working_signal.size());
                y += weights[kk] * working_signal[n - kk];
            }
        }
        SONAR_ALIGNED_REMAINDER_LOOP(k, std::min(M, n + 1), BFP_BLOCK_SIZE) {
            if (k >= M || k > n) continue;
            SONAR_BOUNDS_CHECK(k, weights.size());
            SONAR_BOUNDS_CHECK(n - k, working_signal.size());
            y += weights[k] * working_signal[n - k];
        }
        SONAR_BOUNDS_CHECK(n, equalized.size());
        equalized[n] = y;
    }

    if (use_bfp) {
        double total_scale = 1.0;
        for (const auto& blk : bfp_blocks) {
            total_scale = std::max(total_scale, std::ldexp(1.0, blk.exponent));
        }
        if (total_scale > 1.0) {
            for (auto& s : equalized) {
                s *= total_scale;
            }
        }
    }

    auto t1 = std::chrono::high_resolution_clock::now();

    CMAEqualizerResult result;
    result.equalized_signal = std::move(equalized);
    result.filter_weights = std::move(weights);
    result.iterations_run = convergence_curve.size();
    result.convergence_curve = std::move(convergence_curve);
    result.converged = converged;
    result.final_mse = final_mse;
    result.processing_time_ms =
        std::chrono::duration<double, std::milli>(t1 - t0).count();

    return result;
}

CMAEqualizerResult CMAEqualizer::equalize(const std::vector<double>& real_signal) {
    ComplexVector complex_signal(real_signal.size());
    for (size_t i = 0; i < real_signal.size(); ++i) {
        complex_signal[i] = std::complex<double>(real_signal[i], 0.0);
    }
    return equalize(complex_signal);
}

CMAEqualizerResult CMAEqualizer::equalize_with_doppler(
    const ComplexVector& received_signal,
    const PlatformState& platform,
    const ArrayDeformation& deformation) {
    auto t0 = std::chrono::high_resolution_clock::now();

    if (received_signal.empty()) {
        throw std::invalid_argument("Received signal is empty");
    }

    cached_platform_ = platform;
    cached_deformation_ = deformation;
    resampler_->set_sample_rate(platform.baseband_sample_rate_hz);

    const size_t N = received_signal.size();
    const size_t M = config_.filter_taps;
    if (M == 0) {
        throw std::invalid_argument("filter_taps must be > 0");
    }

    const double R2 = config_.modulus * config_.modulus;
    const double mu = config_.step_size;
    const bool use_bfp = config_.enable_bfp_normalization;
    const bool use_doppler_ff = config_.enable_doppler_feedforward;

    ComplexVector working_signal;

    if (use_doppler_ff) {
        last_doppler_result_ = resampler_->compensate(
            received_signal, platform, deformation, *ray_tracer_, 16);
        working_signal = last_doppler_result_.resampled_signal;
    } else {
        working_signal = received_signal;
    }

    ComplexVector bfp_norm_signal;
    std::vector<BFPBlock> bfp_blocks;

    if (use_bfp) {
        size_t bfp_N = working_signal.size();
        bfp_blocks.reserve((bfp_N + BFP_BLOCK_SIZE - 1) / BFP_BLOCK_SIZE);
        for (size_t bs = 0; bs < bfp_N; bs += BFP_BLOCK_SIZE) {
            const size_t blen = std::min(BFP_BLOCK_SIZE, bfp_N - bs);
            bfp_blocks.push_back(compute_bfp_block(working_signal, bs, blen));
        }
        bfp_norm_signal = reconstruct_bfp_signal(bfp_blocks, bfp_N);
        working_signal = bfp_norm_signal;
    }

    ComplexVector weights(M, std::complex<double>(0.0, 0.0));
    const size_t center_tap = M / 2;
    SONAR_BOUNDS_CHECK(center_tap, weights.size());
    weights[center_tap] = std::complex<double>(1.0, 0.0);

    size_t NN = working_signal.size();
    ComplexVector equalized(NN, std::complex<double>(0.0, 0.0));
    std::vector<double> convergence_curve;
    std::vector<double> doppler_estimates;
    convergence_curve.reserve(config_.max_iterations);
    doppler_estimates.reserve(config_.max_iterations);

    bool converged = false;
    size_t iter = 0;
    double prev_mse = 1e300;
    double final_mse = 0.0;

    ComplexVector signal_buffer(M, std::complex<double>(0.0, 0.0));

    const auto& inst_dop = last_doppler_result_.instantaneous_doppler_hz;
    double fs = std::max(1.0, platform.baseband_sample_rate_hz);

    for (; iter < config_.max_iterations; ++iter) {
        double epoch_mse = 0.0;
        size_t count = 0;
        double epoch_dop = 0.0;

        if (use_bfp && (iter % 10 == 0) && iter > 0) {
            size_t bfp_N = working_signal.size();
            for (size_t bs = 0; bs < bfp_N; bs += BFP_BLOCK_SIZE) {
                const size_t blen = std::min(BFP_BLOCK_SIZE, bfp_N - bs);
                bfp_blocks[bs / BFP_BLOCK_SIZE] =
                    compute_bfp_block(working_signal, bs, blen);
            }
            working_signal = reconstruct_bfp_signal(bfp_blocks, bfp_N);
        }

        std::fill(signal_buffer.begin(), signal_buffer.end(),
                  std::complex<double>(0.0, 0.0));

        for (size_t n = 0; n < NN; ++n) {
            for (size_t k = M - 1; k > 0; --k) {
                SONAR_BOUNDS_CHECK(k, signal_buffer.size());
                SONAR_BOUNDS_CHECK(k - 1, signal_buffer.size());
                signal_buffer[k] = signal_buffer[k - 1];
            }
            SONAR_BOUNDS_CHECK(0, signal_buffer.size());
            SONAR_BOUNDS_CHECK(n, working_signal.size());
            signal_buffer[0] = working_signal[n];

            std::complex<double> y(0.0, 0.0);
            SONAR_ALIGNED_BLOCK_LOOP(k, M, BFP_BLOCK_SIZE) {
                for (size_t j = 0; j < BFP_BLOCK_SIZE; ++j) {
                    const size_t kk = k + j;
                    SONAR_BOUNDS_CHECK(kk, weights.size());
                    SONAR_BOUNDS_CHECK(kk, signal_buffer.size());
                    y += weights[kk] * signal_buffer[kk];
                }
            }
            SONAR_ALIGNED_REMAINDER_LOOP(k, M, BFP_BLOCK_SIZE) {
                SONAR_BOUNDS_CHECK(k, weights.size());
                SONAR_BOUNDS_CHECK(k, signal_buffer.size());
                y += weights[k] * signal_buffer[k];
            }

            SONAR_BOUNDS_CHECK(n, equalized.size());
            equalized[n] = y;

            double mag_sq = std::norm(y);
            if (!std::isfinite(mag_sq) || mag_sq > 1e12) {
                for (size_t kk = 0; kk < M; ++kk) {
                    SONAR_BOUNDS_CHECK(kk, weights.size());
                    weights[kk] = std::complex<double>(0.0, 0.0);
                }
                const size_t center = M / 2;
                SONAR_BOUNDS_CHECK(center, weights.size());
                weights[center] = std::complex<double>(1.0, 0.0);
                std::fill(signal_buffer.begin(), signal_buffer.end(),
                          std::complex<double>(0.0, 0.0));
                continue;
            }
            std::complex<double> e = y * (mag_sq - R2);

            epoch_mse += std::norm(e);
            ++count;

            double dop_bias = 0.0;
            if (use_doppler_ff && !inst_dop.empty()) {
                size_t d_idx = std::min(n, inst_dop.size() - 1);
                dop_bias = inst_dop[d_idx] / fs;
                epoch_dop += inst_dop[d_idx];
            }

            update_weights(weights, e, signal_buffer, mu, dop_bias);
        }

        if (count > 0) {
            epoch_mse /= count;
            epoch_dop /= (double)count;
        }
        convergence_curve.push_back(epoch_mse);
        doppler_estimates.push_back(epoch_dop);
        final_mse = epoch_mse;

        double delta = std::abs(prev_mse - epoch_mse) / (prev_mse + 1e-30);
        if (iter > 10 && delta < config_.convergence_threshold) {
            converged = true;
            break;
        }
        prev_mse = epoch_mse;
    }

    for (size_t n = 0; n < NN; ++n) {
        std::complex<double> y(0.0, 0.0);
        SONAR_ALIGNED_BLOCK_LOOP(k, std::min(M, n + 1), BFP_BLOCK_SIZE) {
            for (size_t j = 0; j < BFP_BLOCK_SIZE; ++j) {
                const size_t kk = k + j;
                if (kk >= M || kk > n) break;
                SONAR_BOUNDS_CHECK(kk, weights.size());
                SONAR_BOUNDS_CHECK(n - kk, working_signal.size());
                y += weights[kk] * working_signal[n - kk];
            }
        }
        SONAR_ALIGNED_REMAINDER_LOOP(k, std::min(M, n + 1), BFP_BLOCK_SIZE) {
            if (k >= M || k > n) continue;
            SONAR_BOUNDS_CHECK(k, weights.size());
            SONAR_BOUNDS_CHECK(n - k, working_signal.size());
            y += weights[k] * working_signal[n - k];
        }
        SONAR_BOUNDS_CHECK(n, equalized.size());
        equalized[n] = y;
    }

    if (use_bfp) {
        double total_scale = 1.0;
        for (const auto& blk : bfp_blocks) {
            total_scale = std::max(total_scale, std::ldexp(1.0, blk.exponent));
        }
        if (total_scale > 1.0) {
            for (auto& s : equalized) {
                s *= total_scale;
            }
        }
    }

    auto t1 = std::chrono::high_resolution_clock::now();

    CMAEqualizerResult result;
    result.equalized_signal = std::move(equalized);
    result.filter_weights = std::move(weights);
    result.iterations_run = convergence_curve.size();
    result.convergence_curve = std::move(convergence_curve);
    result.doppler_estimate_hz = std::move(doppler_estimates);
    result.converged = converged;
    result.final_mse = final_mse;
    result.processing_time_ms =
        std::chrono::duration<double, std::milli>(t1 - t0).count();

    return result;
}

}
