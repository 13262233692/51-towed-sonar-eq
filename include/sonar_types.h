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
#include <cmath>

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
constexpr double SPEED_OF_SOUND_WATER = 1500.0;
constexpr size_t FARROW_POLY_ORDER = 3;

struct Vec3 {
    double x; double y; double z;
    Vec3() : x(0), y(0), z(0) {}
    Vec3(double xx, double yy, double zz) : x(xx), y(yy), z(zz) {}
    inline double norm() const { return std::sqrt(x*x + y*y + z*z); }
    inline double dot(const Vec3& o) const { return x*o.x + y*o.y + z*o.z; }
    inline Vec3 operator+(const Vec3& o) const { return {x+o.x,y+o.y,z+o.z}; }
    inline Vec3 operator-(const Vec3& o) const { return {x-o.x,y-o.y,z-o.z}; }
    inline Vec3 operator*(double s) const { return {x*s,y*s,z*s}; }
};

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
    bool enable_doppler_feedforward = true;
    double doppler_feedback_gain = 0.1;
};

struct CMAEqualizerResult {
    std::vector<std::complex<double>> equalized_signal;
    std::vector<std::complex<double>> filter_weights;
    std::vector<double> convergence_curve;
    std::vector<double> doppler_estimate_hz;
    size_t iterations_run;
    bool converged;
    double final_mse;
    double processing_time_ms;
};

struct SSProfile {
    std::vector<double> depth_m;
    std::vector<double> sound_speed_mps;
};

struct PlatformState {
    double timestamp_s;
    Vec3 sub_position_m;
    Vec3 sub_velocity_mps;
    double sub_yaw_rad;
    double sub_pitch_rad;
    double sub_roll_rad;
    Vec3 target_position_m;
    Vec3 target_velocity_mps;
    double target_rcs_m2;
    double carrier_frequency_hz;
    double baseband_sample_rate_hz;
};

struct ArrayDeformation {
    double timestamp_s;
    std::vector<Vec3> element_positions_m;
    std::vector<Vec3> element_velocities_mps;
    std::array<double, 3> array_bending_curvature;
    std::array<double, 3> array_twist_rate_rps;
};

struct AcousticRay {
    size_t ray_id;
    Vec3 launch_point;
    Vec3 launch_direction;
    double travel_time_s;
    double path_length_m;
    double grazing_angle_rad;
    double incident_angle_rad;
    double reflection_loss_db;
    size_t surface_bounces;
    size_t bottom_bounces;
    double radial_velocity_mps;
    double doppler_shift_hz;
    double complex_weight;
};

struct ResamplingPoly {
    std::array<double, FARROW_POLY_ORDER + 1> coeffs;
    double valid_start_sample;
    double valid_end_sample;
    double rate_ratio;
    double fractional_delay;
};

struct DopplerCompensationResult {
    std::vector<std::complex<double>> resampled_signal;
    std::vector<ResamplingPoly> poly_series;
    std::vector<AcousticRay> active_rays;
    std::vector<double> instantaneous_doppler_hz;
    std::vector<double> radial_velocity_track_mps;
    double total_compensation_gain_db;
    double residual_drift_ppm;
    double processing_time_ms;
};

struct SpectrumArray {
    size_t num_channels;
    size_t fft_size;
    double sample_rate_hz;
    std::vector<double> frequency_axis_hz;
    std::vector<std::vector<double>> power_spectrum_db;
    std::vector<std::vector<std::complex<double>>> complex_spectrum;
    double timestamp_s;
    double integration_time_s;
};

using ComplexVector = std::vector<std::complex<double>>;

}

#endif
