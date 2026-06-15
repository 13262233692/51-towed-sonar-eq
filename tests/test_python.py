import sys
import os
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

try:
    import sonar
except ImportError as e:
    print(f"ERROR: Failed to import sonar module: {e}")
    print(f"Current PYTHONPATH: {os.environ.get('PYTHONPATH', '<not set>')}")
    print(f"Searching in: {sys.path[0]}")
    for f in os.listdir(os.path.dirname(os.path.abspath(__file__))):
        print(f"  Found file: {f}")
    sys.exit(1)

passed = 0
failed = 0

def check(name, condition):
    global passed, failed
    if condition:
        passed += 1
        print(f"  [PASS] {name}")
    else:
        failed += 1
        print(f"  [FAIL] {name}")

def generate_telemetry_bytes(num_frames=1000):
    import struct
    import random
    random.seed(42)

    data = bytearray()
    SYNC = 0xDEADBEEF
    SR = 48000.0
    dt_ns = int(1e9 / SR)

    for f in range(num_frames):
        frame = bytearray()
        frame += struct.pack('<I', SYNC)
        frame += struct.pack('<Q', 1700000000000000000 + f * dt_ns)
        frame += struct.pack('<I', f)
        depth = 100.0 + f * 0.01
        pressure_raw = int(depth / 0.001)
        frame += struct.pack('<H', pressure_raw & 0xFFFF)

        t = f / SR
        for ch in range(sonar.NUM_CHANNELS):
            phase = 2.0 * np.pi * 1000.0 * t + ch * 0.1
            voltage = 1.5 * np.sin(phase) + random.uniform(-0.01, 0.01)
            raw = int(voltage / 3.051850947599719e-05)
            raw = max(-0x800000, min(0x7FFFFF, raw))
            frame += struct.pack('<I', raw & 0xFFFFFF)[:3]

        crc_offset = len(frame)
        crc = 0xFFFF
        for byte in frame:
            crc ^= byte
            for _ in range(8):
                if crc & 1:
                    crc = (crc >> 1) ^ 0xA001
                else:
                    crc >>= 1

        frame += struct.pack('<H', crc & 0xFFFF)
        frame += struct.pack('<B', 0x7E)

        while len(frame) < sonar.PACKET_SIZE:
            frame += b'\x00'
        frame = frame[:sonar.PACKET_SIZE]

        data += frame

    return bytes(data)

def test_module_constants():
    print("\n[Test] Module constants")
    check("NUM_CHANNELS is 64", sonar.NUM_CHANNELS == 64)
    check("PACKET_SIZE > 0", sonar.PACKET_SIZE > 0)
    print(f"  NUM_CHANNELS={sonar.NUM_CHANNELS}, PACKET_SIZE={sonar.PACKET_SIZE}")

def test_hydrophone_unpacker_numpy():
    print("\n[Test] HydrophoneUnpacker (NumPy output)")

    num_frames = 500
    raw_data = generate_telemetry_bytes(num_frames)

    unpacker = sonar.HydrophoneUnpacker()
    check("unpack_to_numpy returns tuple", True)

    result = unpacker.unpack_to_numpy(raw_data)
    timestamps, depths, voltages, valid, corrupt, proc_time = result

    check(f"valid frames ~{num_frames}", valid > num_frames * 0.9)
    check("timestamps shape", timestamps.shape == (valid,))
    check("depths shape", depths.shape == (valid,))
    check("voltages shape", voltages.shape == (valid, sonar.NUM_CHANNELS))
    check("timestamps increasing", np.all(np.diff(timestamps) > 0))
    check("processing time positive", proc_time > 0)

    print(f"  Valid frames: {valid}, Corrupted: {corrupt}")
    print(f"  Processing time: {proc_time:.2f} ms")
    print(f"  Voltage range: [{voltages.min():.4f}, {voltages.max():.4f}] V")

def test_hydrophone_unpacker_threads():
    print("\n[Test] HydrophoneUnpacker thread control")

    unpacker = sonar.HydrophoneUnpacker()
    check("default threads > 0", unpacker.get_num_threads() > 0)

    unpacker.set_num_threads(2)
    check("set_num_threads(2)", unpacker.get_num_threads() == 2)

    est = sonar.HydrophoneUnpacker.estimate_frame_count(sonar.PACKET_SIZE * 100)
    check("estimate_frame_count", est == 100)
    check("num_channels static", unpacker.num_channels() == 64)

def test_cma_equalizer_complex():
    print("\n[Test] CMAEqualizer complex signal")

    np.random.seed(12345)
    N = 2000
    M = 5
    h = np.exp(-0.3 * np.arange(M)) * (0.8 + 0.3 * np.random.randn(M))
    h /= np.sum(np.abs(h))

    t = np.arange(N)
    symbols = np.exp(1j * np.random.uniform(0, 2 * np.pi, N))
    received = np.convolve(symbols, h, mode='full')[:N]
    noise = np.random.randn(N) + 1j * np.random.randn(N)
    snr = 20
    power_sig = np.mean(np.abs(received)**2)
    power_noise = power_sig / (10**(snr / 10))
    received += noise * np.sqrt(power_noise / 2)

    cfg = sonar.CMAConfig()
    cfg.filter_taps = 24
    cfg.step_size = 1e-3
    cfg.modulus = 1.0
    cfg.max_iterations = 100
    cfg.convergence_threshold = 1e-8

    eq = sonar.CMAEqualizer(cfg)
    eq_out, w_out, conv, iters, converged, mse, proc_time = eq.equalize_numpy(received)

    check("equalized length", len(eq_out) == N)
    check("weights length", len(w_out) == cfg.filter_taps)
    check("convergence curve", len(conv) > 0)
    check("iterations in range", 0 < iters <= cfg.max_iterations)
    check("processing time positive", proc_time > 0)
    check("MSE is finite", np.isfinite(mse))

    print(f"  Iterations: {iters}, Converged: {converged}")
    print(f"  Final MSE: {mse:.6e}, Time: {proc_time:.2f} ms")
    if len(conv) >= 2:
        print(f"  First MSE: {conv[0]:.4e}, Last MSE: {conv[-1]:.4e}")

def test_cma_equalizer_real():
    print("\n[Test] CMAEqualizer real signal via standard API")

    N = 1000
    t = np.arange(N)
    signal = np.cos(2 * np.pi * 0.05 * t)
    h = np.array([0.3, 0.7, 0.2, -0.1])
    received = np.convolve(signal, h, mode='full')[:N]
    received += 0.05 * np.random.randn(N)

    eq = sonar.CMAEqualizer()
    cfg = eq.get_config()
    cfg.filter_taps = 16
    cfg.max_iterations = 50
    eq.set_config(cfg)

    received_list = received.tolist()
    result = eq.equalize(received_list)

    check("result type", isinstance(result, sonar.CMAEqualizerResult))
    check("equalized length", len(result.equalized_signal) == N)
    check("weights length", len(result.filter_weights) == cfg.filter_taps)
    check("iterations > 0", result.iterations_run > 0)

    print(f"  Iterations: {result.iterations_run}, Converged: {result.converged}")
    print(f"  Final MSE: {result.final_mse:.6e}")

def test_cma_apply_filter():
    print("\n[Test] CMAEqualizer apply_filter")

    eq = sonar.CMAEqualizer()

    signal = [complex(1.0, 0.0)] * 10
    weights = [complex(0.25, 0.0)] * 4
    result = eq.apply_filter(signal, weights)

    check("output length", len(result) == 10)
    check("first value correct", abs(result[0] - complex(0.25, 0.0)) < 1e-10)
    check("steady-state correct", abs(result[3] - complex(1.0, 0.0)) < 1e-10)

def test_config_classes():
    print("\n[Test] Config classes")

    cfg = sonar.CMAConfig()
    check("default filter_taps", cfg.filter_taps == 32)
    check("default step_size", abs(cfg.step_size - 1e-3) < 1e-12)
    check("default modulus", abs(cfg.modulus - 1.0) < 1e-12)

    cfg.filter_taps = 64
    cfg.step_size = 0.01
    check("filter_taps assignable", cfg.filter_taps == 64)
    check("step_size assignable", abs(cfg.step_size - 0.01) < 1e-12)

    frame = sonar.HydrophoneFrame()
    frame.timestamp_ns = 1234567890
    frame.pressure_depth_m = 150.5
    voltages = list(frame.acoustic_voltage)
    voltages[0] = 1.23
    frame.acoustic_voltage = voltages
    check("frame timestamp", frame.timestamp_ns == 1234567890)
    check("frame depth", abs(frame.pressure_depth_m - 150.5) < 1e-10)
    check("frame voltage", abs(frame.acoustic_voltage[0] - 1.23) < 1e-6)
    check("all channels present", len(frame.acoustic_voltage) == 64)

def main():
    print("=" * 60)
    print("  Python Integration Tests for sonar module")
    print("=" * 60)

    test_module_constants()
    test_hydrophone_unpacker_numpy()
    test_hydrophone_unpacker_threads()
    test_cma_equalizer_complex()
    test_cma_equalizer_real()
    test_cma_apply_filter()
    test_config_classes()

    print("\n" + "=" * 60)
    print(f"  Results: {passed} passed, {failed} failed")
    print("=" * 60)

    sys.exit(0 if failed == 0 else 1)

if __name__ == "__main__":
    main()
