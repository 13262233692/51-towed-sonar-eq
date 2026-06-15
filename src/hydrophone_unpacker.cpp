#include "hydrophone_unpacker.h"
#include "sonar_types.h"
#include <omp.h>
#include <cstring>
#include <chrono>
#include <stdexcept>
#include <stdexcept>

namespace sonar {

namespace {

uint16_t compute_crc16_safe(const uint8_t* data, size_t len) {
    if (data == nullptr && len > 0) {
        throw std::invalid_argument("compute_crc16_safe: null data with nonzero length");
    }
    uint16_t crc = 0xFFFF;
    for (size_t i = 0; i < len; ++i) {
        SONAR_BOUNDS_CHECK(i, len);
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

bool check_sync_word_safe(const uint8_t* ptr, size_t remaining) {
    if (remaining < sizeof(uint32_t)) return false;
    uint32_t sync = 0;
    std::memcpy(&sync, ptr, sizeof(uint32_t));
    return sync == SYNC_WORD;
}

int32_t int24_to_int32_safe(const uint8_t* bytes, size_t remaining) {
    if (remaining < 3) {
        throw std::out_of_range("int24_to_int32_safe: insufficient bytes");
    }
    int32_t val = static_cast<int32_t>(bytes[0])
                | (static_cast<int32_t>(bytes[1]) << 8)
                | (static_cast<int32_t>(bytes[2]) << 16);
    if (val & 0x00800000) {
        val |= 0xFF000000;
    }
    return val;
}

}

HydrophoneUnpacker::HydrophoneUnpacker()
    : num_threads_(0) {
}

void HydrophoneUnpacker::set_num_threads(int num_threads) {
    num_threads_ = num_threads > 0 ? num_threads : omp_get_max_threads();
}

int HydrophoneUnpacker::get_num_threads() const {
    return num_threads_ > 0 ? num_threads_ : omp_get_max_threads();
}

size_t HydrophoneUnpacker::estimate_frame_count(size_t data_size) {
    return data_size / PACKET_SIZE;
}

bool HydrophoneUnpacker::verify_crc16(const RawTelemetryPacket& packet) const {
    const size_t crc_offset = offsetof(RawTelemetryPacket, crc16);
    if (crc_offset >= PACKET_SIZE) {
        throw std::logic_error("verify_crc16: crc_offset out of packet bounds");
    }
    uint16_t computed = compute_crc16_safe(
        reinterpret_cast<const uint8_t*>(&packet), crc_offset);
    return computed == packet.crc16;
}

double HydrophoneUnpacker::convert_pressure(uint16_t raw) const {
    return static_cast<double>(raw) * PRESSURE_SCALE + PRESSURE_OFFSET;
}

double HydrophoneUnpacker::convert_voltage(const uint8_t* raw_bytes) const {
    if (raw_bytes == nullptr) {
        throw std::invalid_argument("convert_voltage: null byte pointer");
    }
    int32_t signed_val = int24_to_int32_safe(raw_bytes, 3);
    return static_cast<double>(signed_val) * VOLTAGE_SCALE;
}

UnpackResult HydrophoneUnpacker::process_chunk(
    const uint8_t* data, size_t start, size_t end) const {

    UnpackResult result;
    result.total_bytes_processed = 0;
    result.valid_frames = 0;
    result.corrupted_frames = 0;
    result.processing_time_ms = 0.0;

    if (data == nullptr) {
        throw std::invalid_argument("process_chunk: null data pointer");
    }
    if (start > end) {
        throw std::invalid_argument("process_chunk: start > end");
    }

    for (size_t pos = start; pos + PACKET_SIZE <= end; pos += PACKET_SIZE) {
        const size_t remaining = end - pos;
        SONAR_BOUNDS_CHECK_INCL(pos + PACKET_SIZE, end);

        if (!check_sync_word_safe(data + pos, remaining)) {
            continue;
        }

        const RawTelemetryPacket* packet =
            reinterpret_cast<const RawTelemetryPacket*>(data + pos);

        if (!verify_crc16(*packet)) {
            result.corrupted_frames++;
            continue;
        }

        HydrophoneFrame frame;
        frame.timestamp_ns = packet->timestamp_ns;
        frame.pressure_depth_m = convert_pressure(packet->pressure_raw);

        for (size_t ch = 0; ch < NUM_HYDROPHONE_CHANNELS; ++ch) {
            SONAR_BOUNDS_CHECK(ch, frame.acoustic_voltage.size());
            const size_t byte_offset = ch * 3;
            SONAR_BOUNDS_CHECK_INCL(byte_offset + 2,
                                    sizeof(packet->channel_data) - 1);
            const uint8_t* ch_ptr = packet->channel_data + byte_offset;
            frame.acoustic_voltage[ch] = convert_voltage(ch_ptr);
        }

        result.frames.push_back(std::move(frame));
        result.valid_frames++;
        result.total_bytes_processed += PACKET_SIZE;
    }

    return result;
}

UnpackResult HydrophoneUnpacker::unpack_stream(
    const uint8_t* data, size_t data_size) const {

    auto t0 = std::chrono::high_resolution_clock::now();

    if (data == nullptr && data_size > 0) {
        throw std::invalid_argument("unpack_stream: null data with nonzero size");
    }

    int threads = num_threads_ > 0 ? num_threads_ : omp_get_max_threads();
    if (threads == 0) threads = 1;
    if (threads > 64) threads = 64;

    size_t chunk_size = 0;
    if (data_size > 0 && static_cast<size_t>(threads) <= data_size) {
        chunk_size = data_size / threads;
        chunk_size = (chunk_size / PACKET_SIZE) * PACKET_SIZE;
        if (chunk_size == 0) chunk_size = PACKET_SIZE;
    } else {
        chunk_size = data_size > 0 ? PACKET_SIZE : 0;
    }

    std::vector<UnpackResult> partial_results(threads);
    for (auto& pr : partial_results) {
        pr.total_bytes_processed = 0;
        pr.valid_frames = 0;
        pr.corrupted_frames = 0;
        pr.processing_time_ms = 0.0;
    }

    #pragma omp parallel for num_threads(threads) \
        default(none) \
        shared(data, data_size, chunk_size, threads, partial_results)
    for (int t = 0; t < threads; ++t) {
        size_t start = 0;
        size_t end = 0;

        if (data_size > 0) {
            start = static_cast<size_t>(t) * chunk_size;
            end = (t == threads - 1) ? data_size : static_cast<size_t>(t + 1) * chunk_size;
            if (start > data_size) start = data_size;
            if (end > data_size) end = data_size;
        }

        UnpackResult local = process_chunk(data, start, end);

        SONAR_BOUNDS_CHECK(static_cast<size_t>(t), partial_results.size());
        partial_results[t] = std::move(local);
    }

    UnpackResult final_result;
    final_result.total_bytes_processed = 0;
    final_result.valid_frames = 0;
    final_result.corrupted_frames = 0;
    final_result.processing_time_ms = 0.0;

    for (size_t t = 0; t < partial_results.size(); ++t) {
        SONAR_BOUNDS_CHECK(t, partial_results.size());
        auto& pr = partial_results[t];
        final_result.frames.insert(
            final_result.frames.end(),
            std::make_move_iterator(pr.frames.begin()),
            std::make_move_iterator(pr.frames.end()));
        final_result.total_bytes_processed += pr.total_bytes_processed;
        final_result.valid_frames += pr.valid_frames;
        final_result.corrupted_frames += pr.corrupted_frames;
    }

    auto t1 = std::chrono::high_resolution_clock::now();
    final_result.processing_time_ms =
        std::chrono::duration<double, std::milli>(t1 - t0).count();

    return final_result;
}

UnpackResult HydrophoneUnpacker::unpack_stream(
    const std::vector<uint8_t>& data) const {
    return unpack_stream(data.data(), data.size());
}

}
