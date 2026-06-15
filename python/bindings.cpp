#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/chrono.h>

#include "hydrophone_unpacker.h"
#include "cma_equalizer.h"
#include "sonar_types.h"

namespace py = pybind11;
using namespace sonar;

namespace {

std::tuple<py::array_t<uint64_t>, py::array_t<double>, py::array_t<double>,
           size_t, size_t, double>
unpack_to_numpy(const HydrophoneUnpacker& unpacker, py::bytes data) {
    std::string s = data;
    UnpackResult result = unpacker.unpack_stream(
        reinterpret_cast<const uint8_t*>(s.data()), s.size());

    size_t n = result.frames.size();

    auto timestamps = py::array_t<uint64_t>(n);
    auto depths = py::array_t<double>(n);
    auto voltages = py::array_t<double>({n, NUM_HYDROPHONE_CHANNELS});

    auto ts_buf = timestamps.mutable_unchecked<1>();
    auto dp_buf = depths.mutable_unchecked<1>();
    auto vg_buf = voltages.mutable_unchecked<2>();

    for (size_t i = 0; i < n; ++i) {
        const auto& f = result.frames[i];
        ts_buf(i) = f.timestamp_ns;
        dp_buf(i) = f.pressure_depth_m;
        for (size_t ch = 0; ch < NUM_HYDROPHONE_CHANNELS; ++ch) {
            vg_buf(i, ch) = f.acoustic_voltage[ch];
        }
    }

    return std::make_tuple(
        timestamps, depths, voltages,
        result.valid_frames,
        result.corrupted_frames,
        result.processing_time_ms
    );
}

py::object cma_equalize_numpy(CMAEqualizer& eq, py::array_t<std::complex<double>> signal) {
    auto buf = signal.unchecked<1>();
    size_t n = buf.shape(0);

    ComplexVector cpp_signal(n);
    for (size_t i = 0; i < n; ++i) {
        cpp_signal[i] = buf(i);
    }

    CMAEqualizerResult result = eq.equalize(cpp_signal);

    auto eq_out = py::array_t<std::complex<double>>(result.equalized_signal.size());
    auto w_out = py::array_t<std::complex<double>>(result.filter_weights.size());
    auto conv_out = py::array_t<double>(result.convergence_curve.size());

    auto eq_buf = eq_out.mutable_unchecked<1>();
    for (size_t i = 0; i < result.equalized_signal.size(); ++i) {
        eq_buf(i) = result.equalized_signal[i];
    }

    auto w_buf = w_out.mutable_unchecked<1>();
    for (size_t i = 0; i < result.filter_weights.size(); ++i) {
        w_buf(i) = result.filter_weights[i];
    }

    auto c_buf = conv_out.mutable_unchecked<1>();
    for (size_t i = 0; i < result.convergence_curve.size(); ++i) {
        c_buf(i) = result.convergence_curve[i];
    }

    return py::make_tuple(
        eq_out, w_out, conv_out,
        result.iterations_run,
        result.converged,
        result.final_mse,
        result.processing_time_ms
    );
}

}

PYBIND11_MODULE(sonar, m) {
    m.doc() = "拖曳线列阵声呐时域基带信号处理与水声通信还原内核";

    py::class_<HydrophoneFrame>(m, "HydrophoneFrame")
        .def(py::init<>())
        .def_readwrite("timestamp_ns", &HydrophoneFrame::timestamp_ns)
        .def_readwrite("pressure_depth_m", &HydrophoneFrame::pressure_depth_m)
        .def_readwrite("acoustic_voltage", &HydrophoneFrame::acoustic_voltage);

    py::class_<UnpackResult>(m, "UnpackResult")
        .def(py::init<>())
        .def_readwrite("frames", &UnpackResult::frames)
        .def_readwrite("total_bytes_processed", &UnpackResult::total_bytes_processed)
        .def_readwrite("valid_frames", &UnpackResult::valid_frames)
        .def_readwrite("corrupted_frames", &UnpackResult::corrupted_frames)
        .def_readwrite("processing_time_ms", &UnpackResult::processing_time_ms);

    py::class_<CMAConfig>(m, "CMAConfig")
        .def(py::init<>())
        .def_readwrite("filter_taps", &CMAConfig::filter_taps)
        .def_readwrite("step_size", &CMAConfig::step_size)
        .def_readwrite("modulus", &CMAConfig::modulus)
        .def_readwrite("leakage", &CMAConfig::leakage)
        .def_readwrite("max_iterations", &CMAConfig::max_iterations)
        .def_readwrite("convergence_threshold", &CMAConfig::convergence_threshold);

    py::class_<CMAEqualizerResult>(m, "CMAEqualizerResult")
        .def(py::init<>())
        .def_readwrite("equalized_signal", &CMAEqualizerResult::equalized_signal)
        .def_readwrite("filter_weights", &CMAEqualizerResult::filter_weights)
        .def_readwrite("convergence_curve", &CMAEqualizerResult::convergence_curve)
        .def_readwrite("iterations_run", &CMAEqualizerResult::iterations_run)
        .def_readwrite("converged", &CMAEqualizerResult::converged)
        .def_readwrite("final_mse", &CMAEqualizerResult::final_mse)
        .def_readwrite("processing_time_ms", &CMAEqualizerResult::processing_time_ms);

    py::class_<HydrophoneUnpacker>(m, "HydrophoneUnpacker")
        .def(py::init<>())
        .def("unpack_stream",
             py::overload_cast<const std::vector<uint8_t>&>(
                 &HydrophoneUnpacker::unpack_stream, py::const_))
        .def("unpack_to_numpy", &unpack_to_numpy)
        .def("set_num_threads", &HydrophoneUnpacker::set_num_threads)
        .def("get_num_threads", &HydrophoneUnpacker::get_num_threads)
        .def_static("estimate_frame_count", &HydrophoneUnpacker::estimate_frame_count)
        .def_static("num_channels", []() { return NUM_HYDROPHONE_CHANNELS; });

    py::class_<CMAEqualizer>(m, "CMAEqualizer")
        .def(py::init<const CMAConfig&>(), py::arg("config") = CMAConfig{})
        .def("equalize", py::overload_cast<const ComplexVector&>(
             &CMAEqualizer::equalize))
        .def("equalize", py::overload_cast<const std::vector<double>&>(
             &CMAEqualizer::equalize))
        .def("equalize_numpy", &cma_equalize_numpy)
        .def("apply_filter", &CMAEqualizer::apply_filter)
        .def("set_config", &CMAEqualizer::set_config)
        .def("get_config", &CMAEqualizer::get_config);

    m.attr("NUM_CHANNELS") = NUM_HYDROPHONE_CHANNELS;
    m.attr("PACKET_SIZE") = PACKET_SIZE;
}
