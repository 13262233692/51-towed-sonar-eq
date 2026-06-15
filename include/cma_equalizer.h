#ifndef CMA_EQUALIZER_H
#define CMA_EQUALIZER_H

#include "sonar_types.h"
#include "doppler_resampler.h"
#include "acoustic_ray_tracer.h"
#include <vector>
#include <complex>
#include <memory>

namespace sonar {

class CMAEqualizer {
public:
    explicit CMAEqualizer(const CMAConfig& config = CMAConfig{});
    ~CMAEqualizer() = default;

    CMAEqualizerResult equalize(const ComplexVector& received_signal);
    CMAEqualizerResult equalize(const std::vector<double>& real_signal);

    CMAEqualizerResult equalize_with_doppler(
        const ComplexVector& received_signal,
        const PlatformState& platform,
        const ArrayDeformation& deformation);

    void set_config(const CMAConfig& config);
    const CMAConfig& get_config() const;

    ComplexVector apply_filter(const ComplexVector& signal,
                               const ComplexVector& weights) const;

    void set_platform_state(const PlatformState& state);
    void set_array_deformation(const ArrayDeformation& deform);
    DopplerCompensationResult get_last_doppler_compensation() const;

private:
    CMAConfig config_;
    PlatformState cached_platform_;
    ArrayDeformation cached_deformation_;
    std::unique_ptr<AcousticRayTracer> ray_tracer_;
    std::unique_ptr<DopplerResampler> resampler_;
    DopplerCompensationResult last_doppler_result_;

    std::complex<double> fir_convolve(const ComplexVector& weights,
                                      const ComplexVector& signal,
                                      size_t index) const;

    void update_weights(ComplexVector& weights,
                        std::complex<double> error,
                        const ComplexVector& signal_buffer,
                        double step_size,
                        double doppler_bias = 0.0) const;
};

}

#endif
