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

void test_cma_thermocline_step_transient() {
    std::cout << "[Test] CMA thermocline step transient (no crash)..." << std::endl;

    const size_t N = 8000;
    const size_t step = N / 2;

    auto signal = generate_step_transient_signal(N, step, 1.0, 100.0, 8.0, 7, 1.0);

    double pre_max = 0.0, post_max = 0.0;
    for (size_t i = 0; i < step; ++i) {
        pre_max = std::max(pre_max, std::abs(signal[i]));
    }
    for (size_t i = step; i < N; ++i) {
        post_max = std::max(post_max, std::abs(signal[i]));
    }
    std::cout << "  Pre-step max amp: " << pre_max
              << ", Post-step max amp: " << post_max << std::endl;

    CMAConfig cfg;
    cfg.filter_taps = 29;
    cfg.step_size = 5e-4;
    cfg.modulus = 1.0;
    cfg.max_iterations = 150;
    cfg.enable_bfp_normalization = true;

    CMAEqualizer eq(cfg);

    bool crashed = false;
    try {
        auto result = eq.equalize(signal);
        TEST_ASSERT(result.equalized_signal.size() == N);
        TEST_ASSERT(result.filter_weights.size() == cfg.filter_taps);
        TEST_ASSERT(result.iterations_run > 0);
        TEST_ASSERT(std::isfinite(result.final_mse));

        for (const auto& v : result.equalized_signal) {
            TEST_ASSERT(std::isfinite(v.real()));
            TEST_ASSERT(std::isfinite(v.imag()));
        }
        for (const auto& w : result.filter_weights) {
            TEST_ASSERT(std::isfinite(w.real()));
            TEST_ASSERT(std::isfinite(w.imag()));
        }

        std::cout << "  Step-transient equalized OK: iterations="
                  << result.iterations_run
                  << ", final MSE=" << std::scientific << result.final_mse << std::endl;
    } catch (const std::exception& e) {
        crashed = true;
        std::cerr << "  EXCEPTION caught: " << e.what() << std::endl;
    }
    TEST_ASSERT(!crashed);
}

void test_cma_aligned_remainder_bounds() {
    std::cout << "[Test] CMA non-aligned remainder boundary guard..." << std::endl;

    for (size_t taps : std::initializer_list<size_t>{7, 9, 15, 17, 23, 25, 31, 33}) {
        const size_t N = 500 + (taps % 8);

        std::vector<double> real_signal(N);
        for (size_t i = 0; i < N; ++i) {
            real_signal[i] = std::sin(0.07 * i) + 0.1 * (i % 13) / 13.0;
        }
        if (N > 200) {
            for (size_t i = 200; i < N; ++i) {
                real_signal[i] *= 50.0;
            }
        }

        CMAConfig cfg;
        cfg.filter_taps = taps;
        cfg.step_size = 2e-4;
        cfg.max_iterations = 30;
        cfg.enable_bfp_normalization = true;

        CMAEqualizer eq(cfg);
        bool ok = false;
        try {
            auto result = eq.equalize(real_signal);
            TEST_ASSERT(result.equalized_signal.size() == N);
            TEST_ASSERT(result.filter_weights.size() == taps);
            ok = true;
        } catch (const std::exception& e) {
            std::cerr << "  taps=" << taps << " threw: " << e.what() << std::endl;
        }
        TEST_ASSERT(ok);
    }
    std::cout << "  All non-aligned tap counts (7-33) processed without bounds violation" << std::endl;
}

void test_cma_bfp_normalization_switch() {
    std::cout << "[Test] CMA BFP normalization on/off consistency..." << std::endl;

    const size_t N = 3000;
    auto signal = generate_step_transient_signal(N, N/2, 1.0, 25.0, 10.0, 5, 1.0);

    CMAConfig cfg_on, cfg_off;
    cfg_on.filter_taps = 16;
    cfg_on.step_size = 1e-3;
    cfg_on.max_iterations = 50;
    cfg_on.enable_bfp_normalization = true;
    cfg_off = cfg_on;
    cfg_off.enable_bfp_normalization = false;

    CMAEqualizer eq_on(cfg_on);
    CMAEqualizer eq_off(cfg_off);

    bool on_ok = true, off_ok = true;
    try {
        auto r_on = eq_on.equalize(signal);
        TEST_ASSERT(r_on.equalized_signal.size() == N);
        TEST_ASSERT(std::isfinite(r_on.final_mse));
    } catch (const std::exception& e) {
        on_ok = false;
        std::cerr << "  BFP=ON threw: " << e.what() << std::endl;
    }

    try {
        auto r_off = eq_off.equalize(signal);
        TEST_ASSERT(r_off.equalized_signal.size() == N);
        TEST_ASSERT(std::isfinite(r_off.final_mse));
    } catch (const std::exception& e) {
        off_ok = false;
        std::cerr << "  BFP=OFF threw: " << e.what() << std::endl;
    }

    TEST_ASSERT(on_ok);
    TEST_ASSERT(off_ok);
    std::cout << "  Both BFP=ON and BFP=OFF completed successfully" << std::endl;
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
    test_cma_thermocline_step_transient();
    test_cma_aligned_remainder_bounds();
    test_cma_bfp_normalization_switch();

    std::cout << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << " Results: " << tests_passed << " passed, "
              << tests_failed << " failed" << std::endl;
    std::cout << "========================================" << std::endl;

    return tests_failed > 0 ? 1 : 0;
}
