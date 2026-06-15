#include "hydrophone_unpacker.h"
#include <omp.h>
#include <cstring>
#include <chrono>
#include <stdexcept>

namespace sonar {

namespace {

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

bool check_sync_word(const uint8_t* ptr) {
    uint32_t sync = 0;
    std::memcpy(&sync, ptr, sizeof(uint32_t));
    return sync == SYNC_WORD;
}

int32_t int24_to_int32(const uint8_t* bytes) {
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
    uint16_t computed = compute_crc16(
        reinterpret_cast<const uint8_t*>(&packet), crc_offset);
    return computed == packet.crc16;
}

double HydrophoneUnpacker::convert_pressure(uint16_t raw) const {
    return static_cast<double>(raw) * PRESSURE_SCALE + PRESSURE_OFFSET;
}

double HydrophoneUnpacker::convert_voltage(const uint8_t* raw_bytes) const {
    int32_t signed_val = int24_to_int32(raw_bytes);
    return static_cast<double>(signed_val) * VOLTAGE_SCALE;
}

UnpackResult HydrophoneUnpacker::process_chunk(
    const uint8_t* data, size_t start, size_t end) const {

    UnpackResult result;
    result.total_bytes_processed = 0;
    result.valid_frames = 0;
    result.corrupted_frames = 0;

    for (size_t pos = start; pos + PACKET_SIZE <= end; pos += PACKET_SIZE) {
        if (!check_sync_word(data + pos)) {
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
            const uint8_t* ch_ptr = packet->channel_data + ch * 3;
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

    int threads = num_threads_ > 0 ? num_threads_ : omp_get_max_threads();
    if (threads == 0) threads = 1;

    size_t chunk_size = data_size / threads;
    chunk_size = (chunk_size / PACKET_SIZE) * PACKET_SIZE;
    if (chunk_size == 0) chunk_size = PACKET_SIZE;

    std::vector<UnpackResult> partial_results(threads);
    for (auto& pr : partial_results) {
        pr.total_bytes_processed = 0;
        pr.valid_frames = 0;
        pr.corrupted_frames = 0;
        pr.processing_time_ms = 0.0;
    }

    #pragma omp parallel for num_threads(threads)
    for (int t = 0; t < threads; ++t) {
        size_t start = t * chunk_size;
        size_t end = (t == threads - 1) ? data_size : (t + 1) * chunk_size;
        partial_results[t] = process_chunk(data, start, end);
    }

    UnpackResult final_result;
    final_result.total_bytes_processed = 0;
    final_result.valid_frames = 0;
    final_result.corrupted_frames = 0;
    for (auto& pr : partial_results) {
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
