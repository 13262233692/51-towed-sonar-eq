#define _USE_MATH_DEFINES
#include "acoustic_ray_tracer.h"
#include "doppler_resampler.h"
#include "spectrum_engine.h"
#include "cma_equalizer.h"
#include "test_utils.h"
#include <cassert>
#include <iostream>
#include <cmath>
#include <chrono>

using namespace sonar;

static size_t g_passed = 0;
static size_t g_failed = 0;

static void check_impl(bool cond, const char* what, const char* file, int line) {
    if (cond) { g_passed++; } else {
        g_failed++;
        std::cout << "  [FAIL] " << what << " at " << file << ":" << line << "\n";
    }
}

#define check(cond, what) check_impl((cond), (what), __FILE__, __LINE__)

static void test_vec3() {
    std::cout << "[Test] Vec3 basic operations...\n";
    Vec3 a(1, 2, 3);
    Vec3 b(4, 5, 6);
    check(std::abs(a.norm() - std::sqrt(14.0)) < 1e-9, "Vec3 norm");
    check(std::abs(a.dot(b) - 32.0) < 1e-9, "Vec3 dot");
    Vec3 c = a + b;
    check(std::abs(c.x - 5.0) < 1e-9, "Vec3 add x");
    Vec3 d = a * 2.0;
    check(std::abs(d.z - 6.0) < 1e-9, "Vec3 mul z");
}

static PlatformState make_high_speed_collision_platform(double t_sec) {
    PlatformState s;
    s.timestamp_s = t_sec;
    s.sub_position_m = Vec3(0, 0, 100);
    double sub_speed = 18.0;
    double vx = std::cos(0.3 * t_sec) * sub_speed;
    double vy = std::sin(0.2 * t_sec) * sub_speed * 0.5;
    double vz = -2.0;
    s.sub_velocity_mps = Vec3(vx, vy, vz);
    s.sub_yaw_rad = 0.3 * t_sec;
    s.sub_pitch_rad = 0.02 * t_sec;
    s.sub_roll_rad = 0.01 * std::sin(t_sec * 0.5);
    s.target_position_m = Vec3(3500, 500, 150);
    double tgt_speed = 22.0;
    s.target_velocity_mps = Vec3(-tgt_speed, -3.0, 0.5);
    s.target_rcs_m2 = 150.0;
    s.carrier_frequency_hz = 3000.0;
    s.baseband_sample_rate_hz = 48000.0;
    return s;
}

static ArrayDeformation make_linear_array_deformation(double t_sec, size_t N = 64) {
    ArrayDeformation a;
    a.timestamp_s = t_sec;
    a.element_positions_m.resize(N);
    a.element_velocities_mps.resize(N);
    for (size_t i = 0; i < N; ++i) {
        double along = (double)i * 0.5;
        double sway = std::sin(t_sec * 0.7 + (double)i * 0.1) * 1.2;
        double heave = std::cos(t_sec * 0.3 + (double)i * 0.05) * 0.5;
        a.element_positions_m[i] = Vec3(along, sway, heave);
        a.element_velocities_mps[i] = Vec3(0,
            0.7 * std::cos(t_sec * 0.7 + (double)i * 0.1) * 1.2,
            -0.3 * std::sin(t_sec * 0.3 + (double)i * 0.05) * 0.5);
    }
    a.array_bending_curvature[0] = 0.001;
    a.array_bending_curvature[1] = 0.0005;
    a.array_bending_curvature[2] = 0.0;
    a.array_twist_rate_rps[0] = 0.0;
    a.array_twist_rate_rps[1] = 0.0;
    a.array_twist_rate_rps[2] = 0.01 * std::sin(t_sec);
    return a;
}

static void test_ray_tracer_basic() {
    std::cout << "[Test] AcousticRayTracer basic SSP interpolation...\n";
    SSProfile ssp;
    ssp.depth_m = {0, 100, 500, 1000, 2000, 4000};
    ssp.sound_speed_mps = {1530, 1510, 1490, 1485, 1500, 1520};
    AcousticRayTracer tracer(ssp);
    double c = tracer.interpolate_sound_speed(300);
    check(c >= 1485 && c <= 1500, "SSP linear interpolation range");
    double c_surf = tracer.interpolate_sound_speed(-10);
    check(std::abs(c_surf - 1530) < 1e-9, "SSP surface clamp");
    double c_deep = tracer.interpolate_sound_speed(10000);
    check(std::abs(c_deep - 1520) < 1e-9, "SSP deep clamp");
    std::cout << "  interpolated c(300m) = " << c << " m/s\n";
}

static void test_ray_tracer_high_speed() {
    std::cout << "[Test] AcousticRayTracer high-speed collision scenario...\n";
    AcousticRayTracer tracer;
    auto state = make_high_speed_collision_platform(0.0);
    auto deform = make_linear_array_deformation(0.0, 16);
    auto rays = tracer.trace_rays(state, deform, 4, 2);
    check(!rays.empty(), "Rays generated");
    std::cout << "  generated " << rays.size() << " rays for grid="
              << (2 * 4 + 1) << "x" << (2 * 4 + 1) << "="
              << (2 * 4 + 1) * (2 * 4 + 1) << " requests\n";
    auto radial = tracer.compute_radial_velocities(state, rays);
    auto doppler = tracer.compute_doppler_shifts(state, rays);
    check(radial.size() == rays.size(), "Radial count matches rays");
    check(doppler.size() == rays.size(), "Doppler count matches rays");
    size_t finite_dop = 0;
    double max_dop = 0;
    for (auto d : doppler) {
        if (std::isfinite(d)) { finite_dop++; max_dop = std::max(max_dop, std::abs(d)); }
    }
    std::cout << "  max |doppler| = " << max_dop << " Hz at fc=" << state.carrier_frequency_hz
              << " Hz, finite=" << finite_dop << "/" << doppler.size() << "\n";
    check(finite_dop > 0, "At least some finite doppler values");
}

static void test_resampler_phase_derotate() {
    std::cout << "[Test] DopplerResampler phase derotation...\n";
    const size_t N = 2048;
    ComplexVector sig(N);
    double fs = 48000;
    double f_sig = 1000;
    double dop = 50;
    for (size_t i = 0; i < N; ++i) {
        double t = (double)i / fs;
        double phase = 2 * M_PI * (f_sig + dop) * t;
        sig[i] = std::complex<double>(std::cos(phase), std::sin(phase));
    }
    DopplerResampler resampler;
    resampler.set_sample_rate(fs);
    std::vector<double> dop_vec(N, dop);
    auto corrected = resampler.phase_derotate(sig, dop_vec, fs);
    check(corrected.size() == N, "Phase derotate size matches");
    double max_phase_error = 0;
    for (size_t i = 100; i < N; ++i) {
        double t = (double)i / fs;
        double expected = 2 * M_PI * f_sig * t;
        std::complex<double> ref(std::cos(expected), std::sin(expected));
        double err = std::abs(corrected[i] * std::conj(ref));
        max_phase_error = std::max(max_phase_error, std::abs(std::arg(corrected[i] * std::conj(ref))));
    }
    std::cout << "  max phase error after derot: " << max_phase_error << " rad\n";
    check(max_phase_error < 0.5, "Phase derotation works approximately correct");
}

static void test_resampler_poly_design() {
    std::cout << "[Test] DopplerResampler polynomial design...\n";
    DopplerResampler resampler;
    const size_t N = 4096;
    double fs = 48000;
    std::vector<double> dop(N);
    double fc = 3000;
    for (size_t i = 0; i < N; ++i) {
        double t = (double)i / fs;
        dop[i] = 80 * std::sin(2 * M_PI * 0.5 * t) + 30;
    }
    auto polys = resampler.design_resampling_polynomials(dop, fc, N);
    check(!polys.empty(), "Polynomial series non-empty");
    std::cout << "  designed " << polys.size() << " polynomials\n";
    size_t all_finite = 0;
    for (const auto& p : polys) {
        bool ok = true;
        for (auto c : p.coeffs) if (!std::isfinite(c)) ok = false;
        if (ok) all_finite++;
    }
    check(all_finite == polys.size(), "All polynomial coeffs finite");
}

static void test_spectrum_engine_welch() {
    std::cout << "[Test] SpectrumEngine Welch PSD...\n";
    const size_t N = 8192;
    const size_t fft_size = 512;
    const double fs = 48000;
    SpectrumEngine engine(fft_size, fs, 2);
    auto axis = engine.get_frequency_axis();
    check(axis.size() == fft_size / 2 + 1, "Freq axis size");
    check(std::abs(axis.back() - fs / 2) < 1e-9, "Nyquist correct");
    std::vector<ComplexVector> channels(2);
    channels[0].resize(N);
    channels[1].resize(N);
    double f1 = 5000;
    double f2 = 12000;
    for (size_t i = 0; i < N; ++i) {
        double t = (double)i / fs;
        channels[0][i] = std::complex<double>(
            std::cos(2*M_PI*f1*t) + 0.3*std::cos(2*M_PI*f2*t),
            std::sin(2*M_PI*f1*t) + 0.3*std::sin(2*M_PI*f2*t));
        channels[1][i] = std::complex<double>(
            0.7 * std::cos(2*M_PI*(f1+200)*t) + std::sin(2*M_PI*20000*t), 0);
    }
    auto sa = engine.compute_welch_array(channels, 0.1);
    check(sa.num_channels == 2, "2 channel count");
    check(sa.power_spectrum_db.size() == 2, "2 PSD channels");
    size_t peak1_bin = 0;
    double peak1 = -1e30;
    for (size_t k = 0; k < sa.power_spectrum_db[0].size(); ++k) {
        if (sa.power_spectrum_db[0][k] > peak1) {
            peak1 = sa.power_spectrum_db[0][k];
            peak1_bin = k;
        }
    }
    double peak_freq = (double)peak1_bin * fs / fft_size;
    std::cout << "  ch0 peak freq = " << peak_freq << " Hz, PSD = " << peak1 << " dB\n";
    check(std::abs(peak_freq - f1) < 1500, "Peak near f1=5kHz");
    check(sa.complex_spectrum.size() == 2, "Complex spectrum channels");
}

static void test_cma_doppler_integration() {
    std::cout << "[Test] CMA with Doppler feed-forward integration (no crash)...\n";
    CMAConfig cfg;
    cfg.filter_taps = 16;
    cfg.max_iterations = 5;
    cfg.step_size = 1e-4;
    cfg.enable_doppler_feedforward = true;
    cfg.enable_bfp_normalization = true;
    CMAEqualizer eq(cfg);
    const size_t N = 512;
    const double fs = 48000;
    const double fc = 3000;
    ComplexVector sig(N);
    double dop = 120;
    for (size_t i = 0; i < N; ++i) {
        double t = (double)i / fs;
        double mod = 1.0 + 0.3 * std::sin(2 * M_PI * 8 * t);
        double phase = 2 * M_PI * (fc / fs * (double)i) + 2 * M_PI * dop * t;
        sig[i] = std::complex<double>(mod * std::cos(phase), mod * std::sin(phase));
    }
    PlatformState state = make_high_speed_collision_platform(0.0);
    state.baseband_sample_rate_hz = fs;
    state.carrier_frequency_hz = fc;
    ArrayDeformation deform = make_linear_array_deformation(0.0, 16);
    auto result = eq.equalize_with_doppler(sig, state, deform);
    check(result.equalized_signal.size() > 0, "Doppler CMA produced output");
    check(result.doppler_estimate_hz.size() == result.iterations_run,
          "Doppler estimate vector length matches iter count");
    std::cout << "  iter=" << result.iterations_run
              << ", final MSE=" << result.final_mse
              << ", eq time=" << result.processing_time_ms << " ms\n";
    auto comp = eq.get_last_doppler_compensation();
    check(comp.active_rays.size() > 0, "Ray tracing ran");
    check(comp.instantaneous_doppler_hz.size() > 0, "Instantaneous doppler vector");
    std::cout << "  resampled size=" << comp.resampled_signal.size()
              << ", rays=" << comp.active_rays.size()
              << ", residual drift=" << comp.residual_drift_ppm << " ppm\n";
}

static void test_full_pipeline_stress() {
    std::cout << "[Test] Full pipeline stress (moving target, 4 channels)...\n";
    SpectrumEngine engine(256, 48000, 4);
    AcousticRayTracer tracer;
    DopplerResampler resampler;
    resampler.set_sample_rate(48000);
    const size_t N = 4096;
    const size_t CH = 4;
    std::vector<ComplexVector> channels(CH);
    for (size_t ch = 0; ch < CH; ++ch) channels[ch].resize(N);
    for (double t_sec = 0; t_sec < 1; t_sec += 0.25) {
        auto state = make_high_speed_collision_platform(t_sec);
        auto deform = make_linear_array_deformation(t_sec, 16);
        auto rays = tracer.trace_rays(state, deform, 2, 1);
        double dop_avg = 0;
        size_t cnt = 0;
        for (const auto& r : rays) {
            if (std::isfinite(r.doppler_shift_hz)) { dop_avg += r.doppler_shift_hz; cnt++; }
        }
        if (cnt) dop_avg /= cnt;
        for (size_t ch = 0; ch < CH; ++ch) {
            double phase_drift = 2 * M_PI * dop_avg * t_sec;
            for (size_t i = 0; i < N; ++i) {
                double tt = t_sec + (double)i / 48000;
                double phase = 2 * M_PI * 3000 * tt + phase_drift + (double)ch * 0.5;
                double amp = 1.0 + 0.5 * std::sin(2 * M_PI * 40 * tt);
                channels[ch][i] = std::complex<double>(amp * std::cos(phase), amp * std::sin(phase));
            }
        }
        auto dop_vec = resampler.estimate_instantaneous_doppler(rays, N, 48000);
        check(dop_vec.size() == N, "Instantaneous doppler length");
    }
    auto sa = engine.compute_welch_array(channels);
    check(sa.num_channels == CH, "Spectrum num channels");
    check(sa.power_spectrum_db.size() == CH, "PSD channel count");
    std::cout << "  final PSD bins/ch OK for " << CH << " channels, " << N << " samples\n";
}

int main() {
    std::cout << "========================================\n";
    std::cout << " Doppler & Spectrum Engine Tests\n";
    std::cout << "========================================\n\n";
    test_vec3();
    test_ray_tracer_basic();
    test_ray_tracer_high_speed();
    test_resampler_phase_derotate();
    test_resampler_poly_design();
    test_spectrum_engine_welch();
    test_cma_doppler_integration();
    test_full_pipeline_stress();
    std::cout << "\n========================================\n";
    std::cout << " Results: " << g_passed << " passed, " << g_failed << " failed\n";
    std::cout << "========================================\n";
    return g_failed > 0 ? 1 : 0;
}
