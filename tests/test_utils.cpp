#define _USE_MATH_DEFINES
#include "test_utils.h"
#include "hydrophone_unpacker.h"
#include <cmath>
#include <random>
#include <cstring>
#include <algorithm>

namespace sonar_test {

uint16_t compute_crc16(const uint8_t* data, size_t len) {
    uint16_t crc = 0xFFFF;
    for (size_t i = 0; i < len; ++i) {
        crc ^= static_cast<uint16_t>(data[i]);
        for (int j = 0; j < 8; ++j) {
            if (crc & 0x0001) {
                crc = (crc >> 1) ^ 0xA001;
            } else {
                crc >>= 1;
            }
        }
    }
    return crc;
}

static void int32_to_int24(int32_t val, uint8_t* out) {
    val = std::clamp(val, -sonar::INT24_MAX - 1, sonar::INT24_MAX);
    out[0] = static_cast<uint8_t>(val & 0xFF);
    out[1] = static_cast<uint8_t>((val >> 8) & 0xFF);
    out[2] = static_cast<uint8_t>((val >> 16) & 0xFF);
}

std::vector<uint8_t> generate_test_telemetry_stream(
    size_t num_frames,
    double signal_freq_hz,
    double sample_rate_hz,
    uint64_t start_timestamp_ns,
    double start_depth_m,
    double depth_increment_m,
    bool insert_corrupted_frames) {

    std::mt19937 rng(42);
    std::uniform_real_distribution<double> noise(-0.01, 0.01);

    const uint64_t dt_ns = static_cast<uint64_t>(1e9 / sample_rate_hz);

    std::vector<uint8_t> stream;
    stream.reserve(num_frames * sonar::PACKET_SIZE);

    for (size_t f = 0; f < num_frames; ++f) {
        sonar::RawTelemetryPacket pkt;
        std::memset(&pkt, 0, sizeof(pkt));

        uint32_t sync = sonar::SYNC_WORD;
        std::memcpy(pkt.sync_word, &sync, sizeof(uint32_t));

        pkt.timestamp_ns = start_timestamp_ns + f * dt_ns;
        pkt.packet_sequence = static_cast<uint32_t>(f);

        double depth = start_depth_m + f * depth_increment_m;
        pkt.pressure_raw = static_cast<uint16_t>(depth / sonar::PRESSURE_SCALE);

        double t = static_cast<double>(f) / sample_rate_hz;
        for (size_t ch = 0; ch < sonar::NUM_HYDROPHONE_CHANNELS; ++ch) {
            double phase = 2.0 * M_PI * signal_freq_hz * t
                         + static_cast<double>(ch) * 0.1;
            double voltage = 1.5 * std::sin(phase) + noise(rng);
            int32_t raw = static_cast<int32_t>(voltage / sonar::VOLTAGE_SCALE);
            int32_to_int24(raw, pkt.channel_data + ch * 3);
        }

        const size_t crc_offset = offsetof(sonar::RawTelemetryPacket, crc16);
        pkt.crc16 = compute_crc16(
            reinterpret_cast<const uint8_t*>(&pkt), crc_offset);

        pkt.frame_delimiter = 0x7E;

        bool corrupt = insert_corrupted_frames && (f % 137 == 0);
        if (corrupt) {
            pkt.crc16 ^= 0xFFFF;
        }

        const uint8_t* bytes = reinterpret_cast<const uint8_t*>(&pkt);
        stream.insert(stream.end(), bytes, bytes + sonar::PACKET_SIZE);
    }

    return stream;
}

std::vector<double> generate_cma_test_signal(
    size_t num_samples,
    double snr_db,
    size_t multipath_taps,
    double carrier_freq) {

    std::mt19937 rng(12345);
    std::normal_distribution<double> gauss(0.0, 1.0);

    std::vector<double> transmit(num_samples);
    for (size_t i = 0; i < num_samples; ++i) {
        double phase = 2.0 * M_PI * carrier_freq * i;
        transmit[i] = std::cos(phase);
    }

    std::vector<double> channel(multipath_taps, 0.0);
    double total = 0.0;
    for (size_t k = 0; k < multipath_taps; ++k) {
        double w = std::exp(-0.3 * k) * (1.0 + 0.3 * gauss(rng));
        channel[k] = w;
        total += std::abs(w);
    }
    for (auto& c : channel) c /= total;

    std::vector<double> received(num_samples, 0.0);
    for (size_t n = 0; n < num_samples; ++n) {
        for (size_t k = 0; k < multipath_taps && k <= n; ++k) {
            received[n] += channel[k] * transmit[n - k];
        }
    }

    double signal_power = 0.0;
    for (auto r : received) signal_power += r * r;
    signal_power /= num_samples;

    double noise_power = signal_power / std::pow(10.0, snr_db / 10.0);
    double noise_std = std::sqrt(noise_power);

    for (auto& r : received) {
        r += noise_std * gauss(rng);
    }

    return received;
}

std::vector<std::complex<double>> generate_complex_cma_signal(
    size_t num_samples,
    double snr_db,
    double modulus) {

    std::mt19937 rng(67890);
    std::normal_distribution<double> gauss(0.0, 1.0);
    std::uniform_real_distribution<double> phase_dist(0.0, 2.0 * M_PI);

    std::vector<std::complex<double>> transmit(num_samples);
    for (size_t i = 0; i < num_samples; ++i) {
        double phi = phase_dist(rng);
        transmit[i] = std::polar(modulus, phi);
    }

    size_t multipath = 5;
    std::vector<std::complex<double>> channel(multipath);
    double total = 0.0;
    for (size_t k = 0; k < multipath; ++k) {
        double mag = std::exp(-0.2 * k) * (0.8 + 0.4 * std::abs(gauss(rng)));
        double ph = 0.5 * k + gauss(rng) * 0.1;
        channel[k] = std::polar(mag, ph);
        total += mag;
    }
    for (auto& c : channel) c /= total;

    std::vector<std::complex<double>> received(num_samples, {0.0, 0.0});
    for (size_t n = 0; n < num_samples; ++n) {
        for (size_t k = 0; k < multipath && k <= n; ++k) {
            received[n] += channel[k] * transmit[n - k];
        }
    }

    double signal_power = 0.0;
    for (const auto& r : received) signal_power += std::norm(r);
    signal_power /= num_samples;

    double noise_power = signal_power / std::pow(10.0, snr_db / 10.0);
    double noise_std = std::sqrt(noise_power / 2.0);

    for (auto& r : received) {
        r += std::complex<double>(noise_std * gauss(rng), noise_std * gauss(rng));
    }

    return received;
}

double compute_mse(const std::vector<double>& a, const std::vector<double>& b) {
    size_t n = std::min(a.size(), b.size());
    if (n == 0) return 0.0;
    double mse = 0.0;
    for (size_t i = 0; i < n; ++i) {
        double d = a[i] - b[i];
        mse += d * d;
    }
    return mse / n;
}

double compute_mse_complex(const std::vector<std::complex<double>>& a,
                            const std::vector<std::complex<double>>& b) {
    size_t n = std::min(a.size(), b.size());
    if (n == 0) return 0.0;
    double mse = 0.0;
    for (size_t i = 0; i < n; ++i) {
        auto d = a[i] - b[i];
        mse += std::norm(d);
    }
    return mse / n;
}

bool approximately_equal(double a, double b, double tol) {
    return std::abs(a - b) < tol;
}

}
