#ifndef TEST_UTILS_H
#define TEST_UTILS_H

#include "sonar_types.h"
#include <vector>
#include <cstdint>
#include <complex>

namespace sonar_test {

uint16_t compute_crc16(const uint8_t* data, size_t len);

std::vector<uint8_t> generate_test_telemetry_stream(
    size_t num_frames,
    double signal_freq_hz = 1000.0,
    double sample_rate_hz = 48000.0,
    uint64_t start_timestamp_ns = 1000000000000ULL,
    double start_depth_m = 100.0,
    double depth_increment_m = 0.01,
    bool insert_corrupted_frames = true
);

std::vector<double> generate_cma_test_signal(
    size_t num_samples,
    double snr_db = 15.0,
    size_t multipath_taps = 5,
    double carrier_freq = 0.05
);

std::vector<std::complex<double>> generate_complex_cma_signal(
    size_t num_samples,
    double snr_db = 15.0,
    double modulus = 1.0
);

double compute_mse(const std::vector<double>& a, const std::vector<double>& b);
double compute_mse_complex(const std::vector<std::complex<double>>& a,
                            const std::vector<std::complex<double>>& b);

bool approximately_equal(double a, double b, double tol = 1e-6);

}

#endif
