#ifndef SONAR_TYPES_H
#define SONAR_TYPES_H

#include <cstdint>
#include <vector>
#include <complex>
#include <array>
#include <string>
#include <chrono>

namespace sonar {

constexpr size_t NUM_HYDROPHONE_CHANNELS = 64;

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
