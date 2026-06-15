#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/chrono.h>

#include "hydrophone_unpacker.h"
#include "cma_equalizer.h"
#include "acoustic_ray_tracer.h"
#include "doppler_resampler.h"
#include "spectrum_engine.h"
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
    auto dop_out = py::array_t<double>(result.doppler_estimate_hz.size());

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

    auto d_buf = dop_out.mutable_unchecked<1>();
    for (size_t i = 0; i < result.doppler_estimate_hz.size(); ++i) {
        d_buf(i) = result.doppler_estimate_hz[i];
    }

    return py::make_tuple(
        eq_out, w_out, conv_out, dop_out,
        result.iterations_run,
        result.converged,
        result.final_mse,
        result.processing_time_ms
    );
}

}

PYBIND11_MODULE(sonar, m) {
    m.doc() = "拖曳线列阵声呐时域基带信号处理与水声通信还原内核";

    py::class_<Vec3>(m, "Vec3")
        .def(py::init<>())
        .def(py::init<double, double, double>())
        .def_readwrite("x", &Vec3::x)
        .def_readwrite("y", &Vec3::y)
        .def_readwrite("z", &Vec3::z)
        .def("norm", &Vec3::norm)
        .def("dot", &Vec3::dot);

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
        .def_readwrite("convergence_threshold", &CMAConfig::convergence_threshold)
        .def_readwrite("enable_bfp_normalization", &CMAConfig::enable_bfp_normalization)
        .def_readwrite("enable_doppler_feedforward", &CMAConfig::enable_doppler_feedforward)
        .def_readwrite("doppler_feedback_gain", &CMAConfig::doppler_feedback_gain);

    py::class_<CMAEqualizerResult>(m, "CMAEqualizerResult")
        .def(py::init<>())
        .def_readwrite("equalized_signal", &CMAEqualizerResult::equalized_signal)
        .def_readwrite("filter_weights", &CMAEqualizerResult::filter_weights)
        .def_readwrite("convergence_curve", &CMAEqualizerResult::convergence_curve)
        .def_readwrite("doppler_estimate_hz", &CMAEqualizerResult::doppler_estimate_hz)
        .def_readwrite("iterations_run", &CMAEqualizerResult::iterations_run)
        .def_readwrite("converged", &CMAEqualizerResult::converged)
        .def_readwrite("final_mse", &CMAEqualizerResult::final_mse)
        .def_readwrite("processing_time_ms", &CMAEqualizerResult::processing_time_ms);

    py::class_<SSProfile>(m, "SSProfile")
        .def(py::init<>())
        .def_readwrite("depth_m", &SSProfile::depth_m)
        .def_readwrite("sound_speed_mps", &SSProfile::sound_speed_mps);

    py::class_<PlatformState>(m, "PlatformState")
        .def(py::init<>())
        .def_readwrite("timestamp_s", &PlatformState::timestamp_s)
        .def_readwrite("sub_position_m", &PlatformState::sub_position_m)
        .def_readwrite("sub_velocity_mps", &PlatformState::sub_velocity_mps)
        .def_readwrite("sub_yaw_rad", &PlatformState::sub_yaw_rad)
        .def_readwrite("sub_pitch_rad", &PlatformState::sub_pitch_rad)
        .def_readwrite("sub_roll_rad", &PlatformState::sub_roll_rad)
        .def_readwrite("target_position_m", &PlatformState::target_position_m)
        .def_readwrite("target_velocity_mps", &PlatformState::target_velocity_mps)
        .def_readwrite("target_rcs_m2", &PlatformState::target_rcs_m2)
        .def_readwrite("carrier_frequency_hz", &PlatformState::carrier_frequency_hz)
        .def_readwrite("baseband_sample_rate_hz", &PlatformState::baseband_sample_rate_hz);

    py::class_<ArrayDeformation>(m, "ArrayDeformation")
        .def(py::init<>())
        .def_readwrite("timestamp_s", &ArrayDeformation::timestamp_s)
        .def_readwrite("element_positions_m", &ArrayDeformation::element_positions_m)
        .def_readwrite("element_velocities_mps", &ArrayDeformation::element_velocities_mps)
        .def_readwrite("array_bending_curvature", &ArrayDeformation::array_bending_curvature)
        .def_readwrite("array_twist_rate_rps", &ArrayDeformation::array_twist_rate_rps);

    py::class_<AcousticRay>(m, "AcousticRay")
        .def(py::init<>())
        .def_readwrite("ray_id", &AcousticRay::ray_id)
        .def_readwrite("launch_point", &AcousticRay::launch_point)
        .def_readwrite("launch_direction", &AcousticRay::launch_direction)
        .def_readwrite("travel_time_s", &AcousticRay::travel_time_s)
        .def_readwrite("path_length_m", &AcousticRay::path_length_m)
        .def_readwrite("grazing_angle_rad", &AcousticRay::grazing_angle_rad)
        .def_readwrite("incident_angle_rad", &AcousticRay::incident_angle_rad)
        .def_readwrite("reflection_loss_db", &AcousticRay::reflection_loss_db)
        .def_readwrite("surface_bounces", &AcousticRay::surface_bounces)
        .def_readwrite("bottom_bounces", &AcousticRay::bottom_bounces)
        .def_readwrite("radial_velocity_mps", &AcousticRay::radial_velocity_mps)
        .def_readwrite("doppler_shift_hz", &AcousticRay::doppler_shift_hz)
        .def_readwrite("complex_weight", &AcousticRay::complex_weight);

    py::class_<ResamplingPoly>(m, "ResamplingPoly")
        .def(py::init<>())
        .def_readwrite("coeffs", &ResamplingPoly::coeffs)
        .def_readwrite("valid_start_sample", &ResamplingPoly::valid_start_sample)
        .def_readwrite("valid_end_sample", &ResamplingPoly::valid_end_sample)
        .def_readwrite("rate_ratio", &ResamplingPoly::rate_ratio)
        .def_readwrite("fractional_delay", &ResamplingPoly::fractional_delay);

    py::class_<DopplerCompensationResult>(m, "DopplerCompensationResult")
        .def(py::init<>())
        .def_readwrite("resampled_signal", &DopplerCompensationResult::resampled_signal)
        .def_readwrite("poly_series", &DopplerCompensationResult::poly_series)
        .def_readwrite("active_rays", &DopplerCompensationResult::active_rays)
        .def_readwrite("instantaneous_doppler_hz", &DopplerCompensationResult::instantaneous_doppler_hz)
        .def_readwrite("radial_velocity_track_mps", &DopplerCompensationResult::radial_velocity_track_mps)
        .def_readwrite("total_compensation_gain_db", &DopplerCompensationResult::total_compensation_gain_db)
        .def_readwrite("residual_drift_ppm", &DopplerCompensationResult::residual_drift_ppm)
        .def_readwrite("processing_time_ms", &DopplerCompensationResult::processing_time_ms);

    py::class_<SpectrumArray>(m, "SpectrumArray")
        .def(py::init<>())
        .def_readwrite("num_channels", &SpectrumArray::num_channels)
        .def_readwrite("fft_size", &SpectrumArray::fft_size)
        .def_readwrite("sample_rate_hz", &SpectrumArray::sample_rate_hz)
        .def_readwrite("frequency_axis_hz", &SpectrumArray::frequency_axis_hz)
        .def_readwrite("power_spectrum_db", &SpectrumArray::power_spectrum_db)
        .def_readwrite("complex_spectrum", &SpectrumArray::complex_spectrum)
        .def_readwrite("timestamp_s", &SpectrumArray::timestamp_s)
        .def_readwrite("integration_time_s", &SpectrumArray::integration_time_s);

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

    py::class_<AcousticRayTracer>(m, "AcousticRayTracer")
        .def(py::init<>())
        .def(py::init<const SSProfile&>())
        .def("set_sound_speed_profile", &AcousticRayTracer::set_sound_speed_profile)
        .def("set_environment", &AcousticRayTracer::set_environment)
        .def("trace_rays", &AcousticRayTracer::trace_rays,
             py::arg("state"), py::arg("array"),
             py::arg("num_rays_per_quadrant") = 16, py::arg("max_bounces") = 2)
        .def("compute_radial_velocities", &AcousticRayTracer::compute_radial_velocities)
        .def("compute_doppler_shifts", &AcousticRayTracer::compute_doppler_shifts)
        .def("interpolate_sound_speed", &AcousticRayTracer::interpolate_sound_speed);

    py::class_<DopplerResampler>(m, "DopplerResampler")
        .def(py::init<>())
        .def("set_sample_rate", &DopplerResampler::set_sample_rate)
        .def("set_farrow_order", &DopplerResampler::set_farrow_order)
        .def("design_resampling_polynomials", &DopplerResampler::design_resampling_polynomials)
        .def("apply_farrow_resampling", &DopplerResampler::apply_farrow_resampling)
        .def("compensate", &DopplerResampler::compensate)
        .def("estimate_instantaneous_doppler", &DopplerResampler::estimate_instantaneous_doppler)
        .def("phase_derotate", &DopplerResampler::phase_derotate);

    py::class_<SpectrumEngine>(m, "SpectrumEngine")
        .def(py::init<>())
        .def(py::init<size_t, double, size_t>(),
             py::arg("fft_size") = 512, py::arg("sample_rate_hz") = 48000.0,
             py::arg("num_channels") = NUM_HYDROPHONE_CHANNELS)
        .def("set_fft_size", &SpectrumEngine::set_fft_size)
        .def("set_sample_rate", &SpectrumEngine::set_sample_rate)
        .def("set_window_type", &SpectrumEngine::set_window_type)
        .def("set_overlap_ratio", &SpectrumEngine::set_overlap_ratio)
        .def("compute_welch_array", &SpectrumEngine::compute_welch_array)
        .def("compute_single_fft", &SpectrumEngine::compute_single_fft)
        .def("get_frequency_axis", &SpectrumEngine::get_frequency_axis)
        .def_static("make_window", &SpectrumEngine::make_window);

    py::class_<CMAEqualizer>(m, "CMAEqualizer")
        .def(py::init<const CMAConfig&>(), py::arg("config") = CMAConfig{})
        .def("equalize", py::overload_cast<const ComplexVector&>(
             &CMAEqualizer::equalize))
        .def("equalize", py::overload_cast<const std::vector<double>&>(
             &CMAEqualizer::equalize))
        .def("equalize_with_doppler", &CMAEqualizer::equalize_with_doppler)
        .def("equalize_numpy", &cma_equalize_numpy)
        .def("apply_filter", &CMAEqualizer::apply_filter)
        .def("set_config", &CMAEqualizer::set_config)
        .def("get_config", &CMAEqualizer::get_config)
        .def("set_platform_state", &CMAEqualizer::set_platform_state)
        .def("set_array_deformation", &CMAEqualizer::set_array_deformation)
        .def("get_last_doppler_compensation", &CMAEqualizer::get_last_doppler_compensation);

    m.attr("NUM_CHANNELS") = NUM_HYDROPHONE_CHANNELS;
    m.attr("PACKET_SIZE") = PACKET_SIZE;
    m.attr("SPEED_OF_SOUND_WATER") = SPEED_OF_SOUND_WATER;
    m.attr("FARROW_POLY_ORDER") = FARROW_POLY_ORDER;
    m.attr("BFP_BLOCK_SIZE") = BFP_BLOCK_SIZE;
}
