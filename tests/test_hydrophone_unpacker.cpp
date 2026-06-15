#include "hydrophone_unpacker.h"
#include "test_utils.h"
#include <iostream>
#include <cassert>
#include <cmath>

using namespace sonar;
using namespace sonar_test;

static int tests_passed = 0;
static int tests_failed = 0;

#define TEST_ASSERT(cond) do { \
    if (!(cond)) { \
        std::cerr << "  FAIL: " << #cond << " at " << __FILE__ << ":" << __LINE__ << std::endl; \
        tests_failed++; \
    } else { \
        tests_passed++; \
    } \
} while(0)

void test_basic_unpack() {
    std::cout << "[Test] Basic unpack functionality..." << std::endl;

    const size_t num_frames = 100;
    auto stream = generate_test_telemetry_stream(num_frames, 1000.0, 48000.0,
                                                 1000000000000ULL, 100.0, 0.01, false);

    HydrophoneUnpacker unpacker;
    auto result = unpacker.unpack_stream(stream);

    TEST_ASSERT(result.valid_frames == num_frames);
    TEST_ASSERT(result.corrupted_frames == 0);
    TEST_ASSERT(result.frames.size() == num_frames);
    TEST_ASSERT(result.total_bytes_processed == num_frames * PACKET_SIZE);
    TEST_ASSERT(result.processing_time_ms > 0.0);

    std::cout << "  Valid frames: " << result.valid_frames
              << ", Time: " << result.processing_time_ms << " ms" << std::endl;
}

void test_timestamp_extraction() {
    std::cout << "[Test] Timestamp nanosecond extraction..." << std::endl;

    const size_t num_frames = 50;
    const uint64_t start_ts = 1700000000000000000ULL;
    const double sr = 96000.0;
    const uint64_t dt_ns = static_cast<uint64_t>(1e9 / sr);

    auto stream = generate_test_telemetry_stream(num_frames, 500.0, sr,
                                                 start_ts, 200.0, 0.0, false);

    HydrophoneUnpacker unpacker;
    auto result = unpacker.unpack_stream(stream);

    TEST_ASSERT(result.valid_frames == num_frames);
    for (size_t i = 0; i < num_frames; ++i) {
        uint64_t expected = start_ts + i * dt_ns;
        TEST_ASSERT(result.frames[i].timestamp_ns == expected);
    }

    std::cout << "  All timestamps correct, start=" << start_ts
              << ", dt=" << dt_ns << "ns" << std::endl;
}

void test_pressure_depth_conversion() {
    std::cout << "[Test] Pressure to depth conversion..." << std::endl;

    const size_t num_frames = 30;
    const double start_depth = 150.0;
    const double inc = 0.5;

    auto stream = generate_test_telemetry_stream(num_frames, 2000.0, 48000.0,
                                                 1000000000ULL, start_depth, inc, false);

    HydrophoneUnpacker unpacker;
    auto result = unpacker.unpack_stream(stream);

    TEST_ASSERT(result.valid_frames == num_frames);
    for (size_t i = 0; i < num_frames; ++i) {
        double expected = start_depth + i * inc;
        double actual = result.frames[i].pressure_depth_m;
        TEST_ASSERT(std::abs(expected - actual) < 0.01);
    }

    std::cout << "  Depth range: " << result.frames.front().pressure_depth_m
              << " -> " << result.frames.back().pressure_depth_m << " m" << std::endl;
}

void test_64_channel_voltage() {
    std::cout << "[Test] 64-channel acoustic voltage extraction..." << std::endl;

    const size_t num_frames = 10;
    auto stream = generate_test_telemetry_stream(num_frames, 800.0, 48000.0,
                                                 0, 50.0, 0.0, false);

    HydrophoneUnpacker unpacker;
    auto result = unpacker.unpack_stream(stream);

    TEST_ASSERT(result.valid_frames == num_frames);
    for (size_t f = 0; f < num_frames; ++f) {
        for (size_t ch = 0; ch < NUM_HYDROPHONE_CHANNELS; ++ch) {
            double v = result.frames[f].acoustic_voltage[ch];
            TEST_ASSERT(v >= -10.0 && v <= 10.0);
        }
        double sum = 0.0;
        for (size_t ch = 0; ch < NUM_HYDROPHONE_CHANNELS; ++ch) {
            sum += result.frames[f].acoustic_voltage[ch];
        }
        TEST_ASSERT(std::abs(sum) < 100.0);
    }

    std::cout << "  All " << NUM_HYDROPHONE_CHANNELS
              << " channels extracted within voltage range" << std::endl;
}

void test_crc_detection() {
    std::cout << "[Test] CRC corrupted frame detection..." << std::endl;

    const size_t num_frames = 500;
    auto stream = generate_test_telemetry_stream(num_frames, 1000.0, 48000.0,
                                                 0, 100.0, 0.01, true);

    HydrophoneUnpacker unpacker;
    auto result = unpacker.unpack_stream(stream);

    size_t expected_corrupt = 0;
    for (size_t f = 0; f < num_frames; ++f) {
        if (f % 137 == 0) expected_corrupt++;
    }

    TEST_ASSERT(result.corrupted_frames == expected_corrupt);
    TEST_ASSERT(result.valid_frames == num_frames - expected_corrupt);

    std::cout << "  Valid: " << result.valid_frames
              << ", Corrupted: " << result.corrupted_frames << std::endl;
}

void test_parallel_accuracy() {
    std::cout << "[Test] OpenMP parallel vs serial consistency..." << std::endl;

    const size_t num_frames = 2000;
    auto stream = generate_test_telemetry_stream(num_frames, 1000.0, 48000.0,
                                                 0, 100.0, 0.01, true);

    HydrophoneUnpacker serial_unpacker;
    serial_unpacker.set_num_threads(1);
    auto serial_result = serial_unpacker.unpack_stream(stream);

    HydrophoneUnpacker parallel_unpacker;
    parallel_unpacker.set_num_threads(4);
    auto parallel_result = parallel_unpacker.unpack_stream(stream);

    TEST_ASSERT(serial_result.valid_frames == parallel_result.valid_frames);
    TEST_ASSERT(serial_result.corrupted_frames == parallel_result.corrupted_frames);
    TEST_ASSERT(serial_result.frames.size() == parallel_result.frames.size());

    for (size_t i = 0; i < serial_result.frames.size(); ++i) {
        TEST_ASSERT(serial_result.frames[i].timestamp_ns ==
                   parallel_result.frames[i].timestamp_ns);
        TEST_ASSERT(approximately_equal(
            serial_result.frames[i].pressure_depth_m,
            parallel_result.frames[i].pressure_depth_m));
        for (size_t ch = 0; ch < NUM_HYDROPHONE_CHANNELS; ++ch) {
            TEST_ASSERT(approximately_equal(
                serial_result.frames[i].acoustic_voltage[ch],
                parallel_result.frames[i].acoustic_voltage[ch]));
        }
    }

    std::cout << "  Serial: " << serial_result.processing_time_ms
              << "ms, Parallel(4): " << parallel_result.processing_time_ms
              << "ms" << std::endl;
}

void test_estimate_frame_count() {
    std::cout << "[Test] Frame count estimation..." << std::endl;

    size_t est = HydrophoneUnpacker::estimate_frame_count(PACKET_SIZE * 1000);
    TEST_ASSERT(est == 1000);

    est = HydrophoneUnpacker::estimate_frame_count(PACKET_SIZE * 5 + 100);
    TEST_ASSERT(est == 5);

    std::cout << "  Estimation correct" << std::endl;
}

int main() {
    std::cout << "========================================" << std::endl;
    std::cout << " Hydrophone Unpacker Unit Tests" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "Packet size: " << PACKET_SIZE << " bytes" << std::endl;
    std::cout << "Channels: " << NUM_HYDROPHONE_CHANNELS << std::endl;
    std::cout << std::endl;

    test_basic_unpack();
    test_timestamp_extraction();
    test_pressure_depth_conversion();
    test_64_channel_voltage();
    test_crc_detection();
    test_parallel_accuracy();
    test_estimate_frame_count();

    std::cout << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << " Results: " << tests_passed << " passed, "
              << tests_failed << " failed" << std::endl;
    std::cout << "========================================" << std::endl;

    return tests_failed > 0 ? 1 : 0;
}
