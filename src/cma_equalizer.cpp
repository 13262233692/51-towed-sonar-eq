#include "cma_equalizer.h"
#include <chrono>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <numeric>

namespace sonar {

CMAEqualizer::CMAEqualizer(const CMAConfig& config)
    : config_(config) {
}

void CMAEqualizer::set_config(const CMAConfig& config) {
    config_ = config;
}

const CMAConfig& CMAEqualizer::get_config() const {
    return config_;
}

std::complex<double> CMAEqualizer::fir_convolve(
    const ComplexVector& weights,
    const ComplexVector& signal,
    size_t index) const {

    std::complex<double> sum(0.0, 0.0);
    const size_t M = weights.size();
    for (size_t k = 0; k < M; ++k) {
        if (index >= k) {
            sum += weights[k] * signal[index - k];
        }
    }
    return sum;
}

void CMAEqualizer::update_weights(
    ComplexVector& weights,
    std::complex<double> error,
    const ComplexVector& signal_buffer,
    double step_size) const {

    const size_t M = weights.size();
    const double leakage = config_.leakage;

    for (size_t k = 0; k < M; ++k) {
        std::complex<double> grad = error * std::conj(signal_buffer[k]);
        weights[k] = (1.0 - leakage) * weights[k] - step_size * grad;
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
    const double R2 = config_.modulus * config_.modulus;
    const double mu = config_.step_size;

    ComplexVector weights(M, std::complex<double>(0.0, 0.0));
    weights[M / 2] = std::complex<double>(1.0, 0.0);

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

        for (size_t n = 0; n < N; ++n) {
            for (size_t k = M - 1; k > 0; --k) {
                signal_buffer[k] = signal_buffer[k - 1];
            }
            signal_buffer[0] = received_signal[n];

            std::complex<double> y(0.0, 0.0);
            for (size_t k = 0; k < M; ++k) {
                y += weights[k] * signal_buffer[k];
            }

            equalized[n] = y;

            double mag_sq = std::norm(y);
            std::complex<double> e = y * (mag_sq - R2);

            epoch_mse += std::norm(e);
            ++count;

            for (size_t k = 0; k < M; ++k) {
                std::complex<double> grad = e * std::conj(signal_buffer[k]);
                weights[k] = (1.0 - config_.leakage) * weights[k] - mu * grad;
            }
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
        for (size_t k = 0; k < M; ++k) {
            if (n >= k) {
                y += weights[k] * received_signal[n - k];
            }
        }
        equalized[n] = y;
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

}
