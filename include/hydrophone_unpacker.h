#ifndef HYDROPHONE_UNPACKER_H
#define HYDROPHONE_UNPACKER_H

#include "sonar_types.h"
#include <cstddef>
#include <cstdint>
#include <vector>

namespace sonar {

#pragma pack(push, 1)
struct RawTelemetryPacket {
    uint8_t  sync_word[4];
    uint64_t timestamp_ns;
    uint32_t packet_sequence;
    uint16_t pressure_raw;
    uint8_t  channel_data[NUM_HYDROPHONE_CHANNELS * 3];
    uint16_t crc16;
    uint8_t  frame_delimiter;
};
#pragma pack(pop)

constexpr size_t PACKET_SIZE = sizeof(RawTelemetryPacket);
constexpr uint32_t SYNC_WORD = 0xDEADBEEF;
constexpr double PRESSURE_SCALE = 0.01;
constexpr double PRESSURE_OFFSET = 0.0;
constexpr double VOLTAGE_SCALE = 3.051850947599719e-05;
constexpr int32_t INT24_MAX = 0x7FFFFF;

class HydrophoneUnpacker {
public:
    HydrophoneUnpacker();
    ~HydrophoneUnpacker() = default;

    UnpackResult unpack_stream(const uint8_t* data, size_t data_size) const;
    UnpackResult unpack_stream(const std::vector<uint8_t>& data) const;

    void set_num_threads(int num_threads);
    int get_num_threads() const;

    static size_t estimate_frame_count(size_t data_size);

private:
    int num_threads_;

    bool verify_crc16(const RawTelemetryPacket& packet) const;
    double convert_pressure(uint16_t raw) const;
    double convert_voltage(const uint8_t* raw_bytes) const;
    UnpackResult process_chunk(const uint8_t* data, size_t start, size_t end) const;
};

}

#endif
