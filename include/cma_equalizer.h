#ifndef CMA_EQUALIZER_H
#define CMA_EQUALIZER_H

#include "sonar_types.h"
#include <vector>
#include <complex>

namespace sonar {

class CMAEqualizer {
public:
    explicit CMAEqualizer(const CMAConfig& config = CMAConfig{});
    ~CMAEqualizer() = default;

    CMAEqualizerResult equalize(const ComplexVector& received_signal);
    CMAEqualizerResult equalize(const std::vector<double>& real_signal);

    void set_config(const CMAConfig& config);
    const CMAConfig& get_config() const;

    ComplexVector apply_filter(const ComplexVector& signal,
                               const ComplexVector& weights) const;

private:
    CMAConfig config_;

    std::complex<double> fir_convolve(const ComplexVector& weights,
                                      const ComplexVector& signal,
                                      size_t index) const;

    void update_weights(ComplexVector& weights,
                        std::complex<double> error,
                        const ComplexVector& signal_buffer,
                        double step_size) const;
};

}

#endif
