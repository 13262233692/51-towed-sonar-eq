#define _USE_MATH_DEFINES
#include "spectrum_engine.h"
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <numeric>

namespace sonar {

SpectrumEngine::SpectrumEngine()
    : fft_size_(512), sample_rate_hz_(48000.0),
      num_channels_(NUM_HYDROPHONE_CHANNELS),
      window_type_("hann"), overlap_ratio_(0.5),
      window_power_correction_(1.0) {
    rebuild_window();
}

SpectrumEngine::SpectrumEngine(size_t fft_size, double sample_rate_hz, size_t num_channels)
    : fft_size_(fft_size), sample_rate_hz_(sample_rate_hz),
      num_channels_(num_channels),
      window_type_("hann"), overlap_ratio_(0.5),
      window_power_correction_(1.0) {
    rebuild_window();
}

void SpectrumEngine::set_fft_size(size_t fft_size) {
    fft_size_ = fft_size;
    rebuild_window();
}

void SpectrumEngine::set_sample_rate(double sample_rate_hz) {
    sample_rate_hz_ = sample_rate_hz;
}

void SpectrumEngine::set_window_type(const std::string& type) {
    window_type_ = type;
    rebuild_window();
}

void SpectrumEngine::set_overlap_ratio(double ratio) {
    overlap_ratio_ = std::max(0.0, std::min(0.95, ratio));
}

std::vector<double> SpectrumEngine::make_window(size_t N, const std::string& type) {
    std::vector<double> w(N, 1.0);
    if (N <= 1) return w;
    if (type == "hann") {
        for (size_t n = 0; n < N; ++n) {
            w[n] = 0.5 * (1.0 - std::cos(2.0 * M_PI * (double)n / (double)(N - 1)));
        }
    } else if (type == "hamming") {
        for (size_t n = 0; n < N; ++n) {
            w[n] = 0.54 - 0.46 * std::cos(2.0 * M_PI * (double)n / (double)(N - 1));
        }
    } else if (type == "blackman") {
        for (size_t n = 0; n < N; ++n) {
            double t = 2.0 * M_PI * (double)n / (double)(N - 1);
            w[n] = 0.42 - 0.5 * std::cos(t) + 0.08 * std::cos(2 * t);
        }
    } else if (type == "blackman-harris") {
        for (size_t n = 0; n < N; ++n) {
            double t = 2.0 * M_PI * (double)n / (double)(N - 1);
            w[n] = 0.35875 - 0.48829 * std::cos(t) + 0.14128 * std::cos(2*t) - 0.01168 * std::cos(3*t);
        }
    }
    return w;
}

void SpectrumEngine::rebuild_window() {
    window_ = make_window(fft_size_, window_type_);
    double win_pow = 0.0, win_sq = 0.0;
    for (size_t n = 0; n < window_.size(); ++n) {
        win_pow += window_[n];
        win_sq += window_[n] * window_[n];
    }
    double N = (double)fft_size_;
    if (win_sq > 1e-12) {
        window_power_correction_ = N / win_sq;
    }
    (void)win_pow;
}

std::vector<double> SpectrumEngine::get_frequency_axis() const {
    std::vector<double> axis(fft_size_ / 2 + 1);
    for (size_t k = 0; k < axis.size(); ++k) {
        axis[k] = (double)k * sample_rate_hz_ / (double)fft_size_;
    }
    return axis;
}

std::vector<std::complex<double>> SpectrumEngine::perform_fft(
    const std::vector<std::complex<double>>& in) {
    size_t N = in.size();
    if (N == 0) return {};
    if ((N & (N - 1)) != 0) {
        size_t np2 = 1;
        while (np2 < N) np2 <<= 1;
        auto padded = in;
        padded.resize(np2, std::complex<double>(0, 0));
        return perform_fft(padded);
    }
    auto buf = in;
    for (size_t i = 1, j = 0; i < N; ++i) {
        size_t bit = N >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) std::swap(buf[i], buf[j]);
    }
    for (size_t len = 2; len <= N; len <<= 1) {
        double ang = -2.0 * M_PI / (double)len;
        std::complex<double> wlen(std::cos(ang), std::sin(ang));
        for (size_t i = 0; i < N; i += len) {
            std::complex<double> w(1, 0);
            for (size_t j = 0; j < len / 2; ++j) {
                auto u = buf[i + j];
                auto v = buf[i + j + len/2] * w;
                buf[i + j] = u + v;
                buf[i + j + len/2] = u - v;
                w *= wlen;
            }
        }
    }
    return buf;
}

double SpectrumEngine::magnitude_to_db(double mag, double ref) {
    if (mag <= 0 || !std::isfinite(mag)) return -200.0;
    return 20.0 * std::log10(mag / ref);
}

SpectrumArray SpectrumEngine::compute_single_fft(
    const ComplexVector& signal, size_t channel_idx) {
    SpectrumArray sa;
    sa.num_channels = 1;
    sa.fft_size = fft_size_;
    sa.sample_rate_hz = sample_rate_hz_;
    sa.frequency_axis_hz = get_frequency_axis();
    sa.timestamp_s = 0.0;
    sa.integration_time_s = (double)fft_size_ / sample_rate_hz_;
    sa.power_spectrum_db.resize(1);
    sa.complex_spectrum.resize(1);
    ComplexVector segment(fft_size_, std::complex<double>(0, 0));
    size_t copy_len = std::min(signal.size(), fft_size_);
    for (size_t i = 0; i < copy_len; ++i) {
        SONAR_BOUNDS_CHECK(i, segment.size());
        SONAR_BOUNDS_CHECK(i, window_.size());
        segment[i] = signal[i] * window_[i];
    }
    auto spectrum = perform_fft(segment);
    sa.complex_spectrum[0].resize(fft_size_ / 2 + 1);
    sa.power_spectrum_db[0].resize(fft_size_ / 2 + 1);
    for (size_t k = 0; k <= fft_size_ / 2; ++k) {
        SONAR_BOUNDS_CHECK(k, sa.complex_spectrum[0].size());
        sa.complex_spectrum[0][k] = spectrum[k] / (double)fft_size_;
        double mag = std::abs(sa.complex_spectrum[0][k]) * window_power_correction_;
        sa.power_spectrum_db[0][k] = magnitude_to_db(mag);
    }
    (void)channel_idx;
    return sa;
}

SpectrumArray SpectrumEngine::compute_welch_array(
    const std::vector<ComplexVector>& channel_signals,
    double integration_time_s) {
    SpectrumArray sa;
    sa.num_channels = channel_signals.size();
    sa.fft_size = fft_size_;
    sa.sample_rate_hz = sample_rate_hz_;
    sa.frequency_axis_hz = get_frequency_axis();
    sa.timestamp_s = 0.0;
    sa.integration_time_s = integration_time_s;
    sa.power_spectrum_db.resize(sa.num_channels);
    sa.complex_spectrum.resize(sa.num_channels);
    size_t spec_bins = fft_size_ / 2 + 1;
    for (size_t ch = 0; ch < sa.num_channels; ++ch) {
        sa.power_spectrum_db[ch].assign(spec_bins, -200.0);
        sa.complex_spectrum[ch].assign(spec_bins, std::complex<double>(0, 0));
    }
    size_t num_segments_total = 0;
    for (size_t ch = 0; ch < channel_signals.size(); ++ch) {
        const auto& sig = channel_signals[ch];
        size_t step = (size_t)std::max((size_t)1,
            (size_t)((double)fft_size_ * (1.0 - overlap_ratio_)));
        std::vector<std::vector<std::complex<double>>> accum_spec(
            spec_bins, std::vector<std::complex<double>>(2, std::complex<double>(0, 0)));
        size_t seg_count = 0;
        for (size_t start = 0; start + fft_size_ <= sig.size(); start += step) {
            ComplexVector segment(fft_size_);
            for (size_t i = 0; i < fft_size_; ++i) {
                SONAR_BOUNDS_CHECK(start + i, sig.size());
                SONAR_BOUNDS_CHECK(i, segment.size());
                SONAR_BOUNDS_CHECK(i, window_.size());
                segment[i] = sig[start + i] * window_[i];
            }
            auto spectrum = perform_fft(segment);
            for (size_t k = 0; k < spec_bins; ++k) {
                std::complex<double> val = spectrum[k] / (double)fft_size_;
                SONAR_BOUNDS_CHECK(k, accum_spec.size());
                accum_spec[k][0] += val;
                accum_spec[k][1] += std::complex<double>(std::norm(val), 0.0);
            }
            seg_count++;
        }
        if (seg_count == 0) {
            size_t copy_len = std::min(sig.size(), fft_size_);
            ComplexVector segment(fft_size_, std::complex<double>(0, 0));
            for (size_t i = 0; i < copy_len; ++i) {
                SONAR_BOUNDS_CHECK(i, segment.size());
                SONAR_BOUNDS_CHECK(i, window_.size());
                segment[i] = sig[i] * window_[i];
            }
            auto spectrum = perform_fft(segment);
            for (size_t k = 0; k < spec_bins; ++k) {
                std::complex<double> val = spectrum[k] / (double)fft_size_;
                SONAR_BOUNDS_CHECK(k, sa.complex_spectrum[ch].size());
                sa.complex_spectrum[ch][k] = val;
                double mag = std::abs(val) * window_power_correction_;
                SONAR_BOUNDS_CHECK(k, sa.power_spectrum_db[ch].size());
                sa.power_spectrum_db[ch][k] = magnitude_to_db(mag);
            }
            num_segments_total++;
        } else {
            for (size_t k = 0; k < spec_bins; ++k) {
                SONAR_BOUNDS_CHECK(k, accum_spec.size());
                sa.complex_spectrum[ch][k] = accum_spec[k][0] / (double)seg_count;
                double avg_power = accum_spec[k][1].real() / (double)seg_count;
                double mag = std::sqrt(std::max(0.0, avg_power)) * window_power_correction_;
                sa.power_spectrum_db[ch][k] = magnitude_to_db(mag);
            }
            num_segments_total += seg_count;
        }
    }
    return sa;
}

}
