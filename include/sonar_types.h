#ifndef SONAR_TYPES_H
#define SONAR_TYPES_H

#include <cstdint>
#include <vector>
#include <complex>
#include <array>
#include <string>
#include <chrono>
#include <cassert>
#include <cstddef>

#define SONAR_ALIGNMENT 8
#define SONAR_PTR_ALIGNED(p) ((((uintptr_t)(p)) & (SONAR_ALIGNMENT - 1)) == 0)

#define SONAR_BOUNDS_CHECK(idx, limit) \
    do { \
        if ((idx) >= (limit)) { \
            assert((idx) < (limit) && "Index out of bounds: missing boundary =="); \
            throw std::out_of_range("sonar: boundary guard violation, idx=" \
                + std::to_string(idx) + " limit=" + std::to_string(limit)); \
        } \
    } while (0)

#define SONAR_BOUNDS_CHECK_INCL(idx, limit) \
    do { \
        if ((idx) > (limit)) { \
            assert((idx) <= (limit) && "Index out of bounds: missing boundary <"); \
            throw std::out_of_range("sonar: inclusive boundary guard violation"); \
        } \
    } while (0)

#define SONAR_ALIGNED_BLOCK_LOOP(i, total, block_size) \
    for (size_t (i) = 0; (i) + (block_size) <= (total); (i) += (block_size))

#define SONAR_ALIGNED_REMAINDER_LOOP(i, total, block_size) \
    for (size_t (i) = ((total) / (block_size)) * (block_size); \
         (i) < (total); ++(i))

namespace sonar {

constexpr size_t NUM_HYDROPHONE_CHANNELS = 64;
constexpr size_t BFP_BLOCK_SIZE = 8;

struct BFPBlock {
    std::vector<std::complex<double>> mantissa;
    int exponent;
};

struct HydrophoneFrame {
    uint64_t timestamp_ns;
    double pressure_depth_m;
    std::array<double, NUM_HYDROPHONE_CHANNELS> acoustic_voltage;
};

struct UnpackResult {
    std::vector<HydrophoneFrame> frames;
    size_t total_bytes_processed;
    size_t valid_frames;
    size_t corrupted_frames;
    double processing_time_ms;
};

struct CMAConfig {
    size_t filter_taps = 32;
    double step_size = 1e-3;
    double modulus = 1.0;
    double leakage = 1e-8;
    size_t max_iterations = 10000;
    double convergence_threshold = 1e-8;
    bool enable_bfp_normalization = true;
};

struct CMAEqualizerResult {
    std::vector<std::complex<double>> equalized_signal;
    std::vector<std::complex<double>> filter_weights;
    std::vector<double> convergence_curve;
    size_t iterations_run;
    bool converged;
    double final_mse;
    double processing_time_ms;
};

struct SSProfile {
    std::vector<double> depth_m;
    std::vector<double> sound_speed_mps;
};

using ComplexVector = std::vector<std::complex<double>>;

}

#endif
