#ifndef SPECTRUM_ENGINE_H
#define SPECTRUM_ENGINE_H

#include "sonar_types.h"
#include <vector>
#include <complex>

namespace sonar {

class SpectrumEngine {
public:
    SpectrumEngine();
    SpectrumEngine(size_t fft_size, double sample_rate_hz,
                   size_t num_channels = NUM_HYDROPHONE_CHANNELS);

    void set_fft_size(size_t fft_size);
    void set_sample_rate(double sample_rate_hz);
    void set_window_type(const std::string& type);
    void set_overlap_ratio(double ratio);

    SpectrumArray compute_welch_array(
        const std::vector<ComplexVector>& channel_signals,
        double integration_time_s = 0.5);

    SpectrumArray compute_single_fft(
        const ComplexVector& signal,
        size_t channel_idx = 0);

    std::vector<double> get_frequency_axis() const;

    static std::vector<double> make_window(size_t N, const std::string& type);

private:
    size_t fft_size_;
    double sample_rate_hz_;
    size_t num_channels_;
    std::string window_type_;
    double overlap_ratio_;
    std::vector<double> window_;
    double window_power_correction_;

    void rebuild_window();
    static std::vector<std::complex<double>> perform_fft(
        const std::vector<std::complex<double>>& in);
    static double magnitude_to_db(double mag, double ref = 1.0);
};

}

#endif
