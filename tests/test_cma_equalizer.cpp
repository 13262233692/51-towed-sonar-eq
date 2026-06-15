#include "cma_equalizer.h"
#include "test_utils.h"
#include <iostream>
#include <cassert>
#include <cmath>
#include <algorithm>
#include <numeric>

using namespace sonar;
using namespace sonar_test;

static int tests_passed = 0;
static int tests_failed = 0;

#define TEST_ASSERT(cond) do { \
    if (!(cond)) { \
        std::cerr << "  FAIL: " << #cond << " at " << __FILE__ << ":" << __LINE__ << std::endl; \
        tests_failed++; \
    } else { \
        tests_passed++; \
    } \
} while(0)

void test_cma_basic_real_signal() {
    std::cout << "[Test] CMA basic real signal equalization..." << std::endl;

    const size_t N = 5000;
    auto received = generate_cma_test_signal(N, 20.0, 4, 0.05);

    CMAConfig cfg;
    cfg.filter_taps = 32;
    cfg.step_size = 1e-4;
    cfg.max_iterations = 200;
    cfg.convergence_threshold = 1e-7;

    CMAEqualizer eq(cfg);
    auto result = eq.equalize(received);

    TEST_ASSERT(result.equalized_signal.size() == N);
    TEST_ASSERT(result.filter_weights.size() == cfg.filter_taps);
    TEST_ASSERT(result.iterations_run > 0);
    TEST_ASSERT(result.iterations_run <= cfg.max_iterations);
    TEST_ASSERT(result.processing_time_ms > 0.0);

    double input_power = 0.0;
    for (auto v : received) input_power += v * v;
    input_power /= N;

    double output_power = 0.0;
    for (const auto& v : result.equalized_signal) {
        output_power += std::norm(v);
    }
    output_power /= N;

    TEST_ASSERT(output_power > 1e-6);

    std::cout << "  Iterations: " << result.iterations_run
              << ", Converged: " << (result.converged ? "YES" : "NO")
              << ", Final MSE: " << result.final_mse << std::endl;
    std::cout << "  Input power: " << input_power
              << ", Output power: " << output_power << std::endl;
}

void test_cma_complex_signal() {
    std::cout << "[Test] CMA complex PSK signal equalization..." << std::endl;

    const size_t N = 3000;
    auto received = generate_complex_cma_signal(N, 20.0, 1.0);

    CMAConfig cfg;
    cfg.filter_taps = 24;
    cfg.step_size = 5e-4;
    cfg.modulus = 1.0;
    cfg.max_iterations = 150;
    cfg.convergence_threshold = 1e-8;

    CMAEqualizer eq(cfg);
    auto result = eq.equalize(received);

    TEST_ASSERT(result.equalized_signal.size() == N);
    TEST_ASSERT(result.filter_weights.size() == cfg.filter_taps);
    TEST_ASSERT(result.convergence_curve.size() == result.iterations_run);

    bool error_decreasing = true;
    if (result.convergence_curve.size() >= 3) {
        double avg_first = 0.0;
        double avg_last = 0.0;
        size_t n = std::min(size_t(10), result.convergence_curve.size() / 3);
        for (size_t i = 0; i < n; ++i) {
            avg_first += result.convergence_curve[i];
            avg_last += result.convergence_curve[result.convergence_curve.size() - 1 - i];
        }
        avg_first /= n;
        avg_last /= n;
        error_decreasing = (avg_last < avg_first * 1.5);
    }
    TEST_ASSERT(error_decreasing);

    double mean_mag = 0.0;
    for (const auto& v : result.equalized_signal) {
        mean_mag += std::abs(v);
    }
    mean_mag /= N;
    TEST_ASSERT(mean_mag > 0.1 && mean_mag < 5.0);

    std::cout << "  Iterations: " << result.iterations_run
              << ", Converged: " << (result.converged ? "YES" : "NO")
              << ", Final MSE: " << std::scientific << result.final_mse << std::endl;
}

void test_cma_fir_convolution() {
    std::cout << "[Test] FIR filter convolution..." << std::endl;

    CMAEqualizer eq;

    ComplexVector weights(4, std::complex<double>(0.25, 0.0));

    ComplexVector signal(10, std::complex<double>(0.0, 0.0));
    for (size_t i = 0; i < signal.size(); ++i) {
        signal[i] = std::complex<double>(static_cast<double>(i + 1), 0.0);
    }

    auto filtered = eq.apply_filter(signal, weights);

    TEST_ASSERT(filtered.size() == signal.size());
    TEST_ASSERT(std::abs(filtered[0].real() - 0.25) < 1e-10);
    TEST_ASSERT(std::abs(filtered[1].real() - 0.75) < 1e-10);
    TEST_ASSERT(std::abs(filtered[2].real() - 1.5) < 1e-10);
    TEST_ASSERT(std::abs(filtered[3].real() - 2.5) < 1e-10);
    TEST_ASSERT(std::abs(filtered[4].real() - 3.5) < 1e-10);

    std::cout << "  FIR convolution correct" << std::endl;
}

void test_cma_config_getset() {
    std::cout << "[Test] CMA config getter/setter..." << std::endl;

    CMAConfig cfg;
    cfg.filter_taps = 64;
    cfg.step_size = 0.005;
    cfg.modulus = 2.0;
    cfg.leakage = 1e-6;
    cfg.max_iterations = 5000;
    cfg.convergence_threshold = 1e-10;

    CMAEqualizer eq(cfg);
    auto got_cfg = eq.get_config();

    TEST_ASSERT(got_cfg.filter_taps == 64);
    TEST_ASSERT(approximately_equal(got_cfg.step_size, 0.005, 1e-12));
    TEST_ASSERT(approximately_equal(got_cfg.modulus, 2.0));
    TEST_ASSERT(approximately_equal(got_cfg.leakage, 1e-6));
    TEST_ASSERT(got_cfg.max_iterations == 5000);
    TEST_ASSERT(approximately_equal(got_cfg.convergence_threshold, 1e-10));

    CMAConfig cfg2;
    cfg2.filter_taps = 16;
    eq.set_config(cfg2);
    auto got2 = eq.get_config();
    TEST_ASSERT(got2.filter_taps == 16);

    std::cout << "  Config getter/setter correct" << std::endl;
}

void test_cma_edge_cases() {
    std::cout << "[Test] CMA edge cases..." << std::endl;

    CMAEqualizer eq;

    bool threw = false;
    try {
        eq.equalize(ComplexVector());
    } catch (const std::invalid_argument&) {
        threw = true;
    } catch (...) {
        threw = true;
    }
    TEST_ASSERT(threw);

    const size_t N = 100;
    ComplexVector ones(N, std::complex<double>(1.0, 0.0));

    CMAConfig cfg;
    cfg.filter_taps = 8;
    cfg.step_size = 1e-5;
    cfg.max_iterations = 5;
    eq.set_config(cfg);

    auto result = eq.equalize(ones);
    TEST_ASSERT(result.equalized_signal.size() == N);
    TEST_ASSERT(result.iterations_run == 5);

    std::cout << "  Edge cases handled correctly" << std::endl;
}

void test_cma_convergence_curve() {
    std::cout << "[Test] CMA convergence curve generation..." << std::endl;

    const size_t N = 2000;
    auto received = generate_complex_cma_signal(N, 25.0, 1.0);

    CMAConfig cfg;
    cfg.filter_taps = 16;
    cfg.step_size = 2e-3;
    cfg.max_iterations = 50;
    cfg.convergence_threshold = 1e-9;

    CMAEqualizer eq(cfg);
    auto result = eq.equalize(received);

    TEST_ASSERT(!result.convergence_curve.empty());
    TEST_ASSERT(result.convergence_curve.size() <= cfg.max_iterations);

    for (auto v : result.convergence_curve) {
        TEST_ASSERT(std::isfinite(v));
        TEST_ASSERT(v >= 0.0);
    }

    if (result.convergence_curve.size() >= 2) {
        double first_half_avg = 0.0;
        double second_half_avg = 0.0;
        size_t half = result.convergence_curve.size() / 2;
        for (size_t i = 0; i < half; ++i) {
            first_half_avg += result.convergence_curve[i];
        }
        for (size_t i = half; i < result.convergence_curve.size(); ++i) {
            second_half_avg += result.convergence_curve[i];
        }
        first_half_avg /= half;
        second_half_avg /= (result.convergence_curve.size() - half);
        std::cout << "  Avg MSE first half: " << std::scientific << first_half_avg
                  << ", second half: " << second_half_avg << std::endl;
    }
}

int main() {
    std::cout << "========================================" << std::endl;
    std::cout << " CMA Blind Equalizer Unit Tests" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << std::endl;

    test_cma_basic_real_signal();
    test_cma_complex_signal();
    test_cma_fir_convolution();
    test_cma_config_getset();
    test_cma_edge_cases();
    test_cma_convergence_curve();

    std::cout << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << " Results: " << tests_passed << " passed, "
              << tests_failed << " failed" << std::endl;
    std::cout << "========================================" << std::endl;

    return tests_failed > 0 ? 1 : 0;
}
