#ifndef DOPPLER_RESAMPLER_H
#define DOPPLER_RESAMPLER_H

#include "sonar_types.h"
#include "acoustic_ray_tracer.h"
#include <vector>
#include <complex>

namespace sonar {

class DopplerResampler {
public:
    DopplerResampler();

    void set_sample_rate(double sample_rate_hz);
    void set_farrow_order(size_t order);

    std::vector<ResamplingPoly> design_resampling_polynomials(
        const std::vector<double>& instantaneous_doppler_hz,
        double carrier_frequency_hz,
        size_t signal_length,
        double poly_update_interval_samples = 64.0);

    std::vector<std::complex<double>> apply_farrow_resampling(
        const std::vector<std::complex<double>>& input,
        const std::vector<ResamplingPoly>& polys);

    DopplerCompensationResult compensate(
        const std::vector<std::complex<double>>& input_signal,
        const PlatformState& state,
        const ArrayDeformation& deformation,
        AcousticRayTracer& tracer,
        size_t num_rays = 16);

    std::vector<double> estimate_instantaneous_doppler(
        const std::vector<AcousticRay>& rays,
        size_t signal_length,
        double sample_rate_hz);

    ComplexVector phase_derotate(
        const ComplexVector& signal,
        const std::vector<double>& doppler_hz,
        double sample_rate_hz);

private:
    double sample_rate_hz_;
    size_t farrow_order_;
    static constexpr size_t FIR_HALF_TAPS = 6;

    double lagrange_kernel(double fractional, size_t order, size_t tap);
    std::vector<std::vector<double>> compute_farrow_coeff_table(double rate_ratio);
    std::complex<double> farrow_interpolate(
        const std::vector<std::complex<double>>& buf,
        double mu,
        const std::vector<std::vector<double>>& coeffs);
};

}

#endif
