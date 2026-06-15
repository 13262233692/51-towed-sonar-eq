#include "hydrophone_unpacker.h"
#include "cma_equalizer.h"
#include "test_utils.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>

using namespace sonar;
using namespace sonar_test;

int main(int argc, char* argv[]) {
    std::string output_dir = ".";
    size_t num_frames = 10000;
    size_t num_samples = 10000;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--frames" && i + 1 < argc) {
            num_frames = std::stoull(argv[++i]);
        } else if (arg == "--samples" && i + 1 < argc) {
            num_samples = std::stoull(argv[++i]);
        } else if (arg == "--out" && i + 1 < argc) {
            output_dir = argv[++i];
        }
    }

    std::cout << "Generating test data..." << std::endl;
    std::cout << "  Frames for telemetry: " << num_frames << std::endl;
    std::cout << "  Samples for CMA: " << num_samples << std::endl;
    std::cout << std::endl;

    {
        auto stream = generate_test_telemetry_stream(num_frames, 1000.0, 48000.0,
                                                     1700000000000000000ULL,
                                                     150.0, 0.005, true);
        std::string path = output_dir + "/telemetry_stream.bin";
        std::ofstream ofs(path, std::ios::binary);
        if (!ofs) {
            std::cerr << "ERROR: Cannot open " << path << " for writing" << std::endl;
            return 1;
        }
        ofs.write(reinterpret_cast<const char*>(stream.data()), stream.size());
        std::cout << "Written: " << path << " (" << stream.size() << " bytes)" << std::endl;

        HydrophoneUnpacker unpacker;
        auto result = unpacker.unpack_stream(stream);
        std::cout << "  Unpacked: " << result.valid_frames << " valid, "
                  << result.corrupted_frames << " corrupted" << std::endl;
    }

    {
        auto signal = generate_complex_cma_signal(num_samples, 20.0, 1.0);
        std::string path = output_dir + "/cma_signal.csv";
        std::ofstream ofs(path);
        if (!ofs) {
            std::cerr << "ERROR: Cannot open " << path << " for writing" << std::endl;
            return 1;
        }
        ofs << "index,real,imag,magnitude\n";
        for (size_t i = 0; i < signal.size(); ++i) {
            ofs << i << ","
                << signal[i].real() << ","
                << signal[i].imag() << ","
                << std::abs(signal[i]) << "\n";
        }
        std::cout << "Written: " << path << " (" << signal.size() << " samples)" << std::endl;
    }

    {
        auto signal = generate_cma_test_signal(num_samples, 15.0, 5, 0.03);
        std::string path = output_dir + "/cma_signal_real.csv";
        std::ofstream ofs(path);
        if (!ofs) {
            std::cerr << "ERROR: Cannot open " << path << " for writing" << std::endl;
            return 1;
        }
        ofs << "index,value\n";
        for (size_t i = 0; i < signal.size(); ++i) {
            ofs << i << "," << signal[i] << "\n";
        }
        std::cout << "Written: " << path << " (" << signal.size() << " samples)" << std::endl;
    }

    std::cout << "\nAll test data generated successfully." << std::endl;
    return 0;
}
