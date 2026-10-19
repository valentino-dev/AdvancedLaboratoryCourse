"""Analyze the supplied GaAs transient-grating and spectrum data.

The delay units and the physical meaning of the grating-file suffixes are not
documented in the data files. Results involving these values are therefore
exploratory; see the lab report for the associated limitations.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import linregress
from scipy.ndimage import uniform_filter1d


PROJECT_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = PROJECT_DIR / "Data TG 4_bis"
PLOT_DIR = PROJECT_DIR / "plots"
PLOT_DIR.mkdir(exist_ok=True)

TRACE_FILES = (
    "GaAs_RT_PP_X.dat",
    "GaAs_RT_5.6_X.dat",
    "GaAs_RT_6.5_X.dat",
    "GaAs_RT_8_X.dat",
    "GaAs_RT_12_X.dat",
)

FIT_START = 0.5
FIT_END = 15.0
ASSUMED_PERIOD_UNIT_METRES = 1e-6  # Exploratory assumption: suffix is in micrometres.

plt.style.use("seaborn-v0_8-whitegrid")


def load_trace(filename):
    """Load a trace: first column is delay, remaining columns are repeats."""
    return np.loadtxt(DATA_DIR / filename)


def exponential_with_offset(delay, amplitude, lifetime, offset):
    """Single-exponential decay model used for exploratory fits."""
    return offset + amplitude * np.exp(-delay / lifetime)


def gaussian_with_offset(wavelength, amplitude, center, sigma, offset):
    """Gaussian spectral peak with a constant background."""
    return offset + amplitude * np.exp(-0.5 * ((wavelength - center) / sigma) ** 2)


def interferometric_autocorrelation_model(
    time,
    offset,
    dc_amplitude,
    omega_amplitude,
    omega_phase,
    twice_omega_amplitude,
    twice_omega_phase,
    center_time,
    sigma_time,
    fringe_frequency,
):
    """Gaussian-envelope IAC fit with the tutorial's omega and 2-omega terms."""
    relative_time = time - center_time
    envelope = np.exp(-0.5 * (relative_time / sigma_time) ** 2)
    fundamental = omega_amplitude * np.cos(
        2 * np.pi * fringe_frequency * relative_time + omega_phase
    )
    second_harmonic = twice_omega_amplitude * np.cos(
        4 * np.pi * fringe_frequency * relative_time + twice_omega_phase
    )
    return offset + envelope * (dc_amplitude + fundamental + second_harmonic)


def fit_trace(delay, signal):
    """Fit a single exponential after time zero and return parameters/mask."""
    fit_mask = (delay >= FIT_START) & (delay <= min(FIT_END, delay.max()))
    initial_guess = (
        signal[fit_mask][0] - signal[fit_mask][-1],
        3.0,
        signal[fit_mask][-1],
    )
    parameters, _ = curve_fit(
        exponential_with_offset,
        delay[fit_mask],
        signal[fit_mask],
        p0=initial_guess,
        maxfev=20_000,
    )
    return parameters, fit_mask


def nominal_period_and_wave_number(filename):
    """Read the filename suffix as an assumed period and calculate q=2*pi/L."""
    if "_PP_" in filename:
        return 0.0, 0.0

    period_label = float(filename.split("_")[2])
    period_metres = period_label * ASSUMED_PERIOD_UNIT_METRES
    wave_number = 2 * np.pi / period_metres
    return period_label, wave_number


def analyze_delay_traces():
    """Plot repeated delay traces, fitted curves, and exploratory rate vs q^2."""
    figure, axis = plt.subplots(figsize=(8, 5))
    fit_results = []

    for filename in TRACE_FILES:
        data = load_trace(filename)
        delay = data[:, 0]
        repeats = data[:, 1:]
        mean_signal = repeats.mean(axis=1)
        standard_error = repeats.std(axis=1, ddof=1) / np.sqrt(repeats.shape[1])

        parameters, fit_mask = fit_trace(delay, mean_signal)
        nominal_period, wave_number = nominal_period_and_wave_number(filename)
        lifetime = parameters[1]
        fit_results.append((filename, nominal_period, wave_number, lifetime))

        label = filename.removeprefix("GaAs_RT_").removesuffix("_X.dat")
        axis.errorbar(
            delay,
            mean_signal,
            yerr=standard_error,
            fmt=".",
            markersize=3,
            label=label,
        )
        axis.plot(
            delay[fit_mask],
            exponential_with_offset(delay[fit_mask], *parameters),
            linewidth=1.5,
        )

    axis.set(
        xlabel="Delay (dataset units)",
        ylabel="Signal (a.u.)",
        title="GaAs pump-probe and transient-grating traces",
    )
    axis.legend()
    figure.tight_layout()
    figure.savefig(PLOT_DIR / "gaas_decay.png", dpi=180)
    plt.close(figure)

    wave_numbers = np.array([result[2] for result in fit_results])
    lifetimes = np.array([result[3] for result in fit_results])
    rates = 1 / lifetimes
    slope, intercept, r_value, _, _ = linregress(wave_numbers**2, rates)

    q_squared = wave_numbers**2
    q_squared_in_1e12_per_m2 = q_squared * 1e-12
    figure, axis = plt.subplots(figsize=(6, 4))
    axis.scatter(q_squared_in_1e12_per_m2, rates, label="Fits")
    x_line = np.linspace(0, q_squared_in_1e12_per_m2.max(), 100)
    axis.plot(x_line, intercept + slope * x_line * 1e12, label="Linear fit")
    axis.set(
        xlabel=r"$q^2$ ($10^{12}$ m$^{-2}$)",
        ylabel="Decay rate (dataset units$^{-1}$)",
        title="Rate versus grating wave number squared",
    )
    axis.legend()
    figure.tight_layout()
    figure.savefig(PLOT_DIR / "rate_vs_q2.png", dpi=180)
    plt.close(figure)

    return fit_results, slope, intercept, r_value**2


def analyze_spectrum():
    """Fit and plot the spectrum; return peak, FWHM and baseline information."""
    spectrum = np.loadtxt(DATA_DIR / "spektrum.csv", delimiter=",")
    wavelength = spectrum[:, 0]
    intensity = spectrum[:, 1]

    median_intensity = np.median(intensity)
    peak_index = np.argmax(intensity)
    initial_guess = (
        intensity.max() - median_intensity,
        wavelength[peak_index],
        5.0,
        median_intensity,
    )
    parameters, _ = curve_fit(
        gaussian_with_offset,
        wavelength,
        intensity,
        p0=initial_guess,
        maxfev=20_000,
    )
    amplitude, center, sigma, offset = parameters
    fitted_peak_intensity = offset + amplitude
    half_maximum = offset + amplitude / 2
    fwhm = 2 * np.sqrt(2 * np.log(2)) * abs(sigma)
    left_half_max = center - fwhm / 2
    right_half_max = center + fwhm / 2

    figure, axis = plt.subplots(figsize=(7, 4))
    axis.plot(wavelength, intensity, label="Measured spectrum")
    axis.plot(
        wavelength,
        gaussian_with_offset(wavelength, *parameters),
        label="Gaussian fit",
    )
    axis.axhline(
        median_intensity,
        color="tab:green",
        linestyle="--",
        label=f"Median = {median_intensity:.1f}",
    )
    axis.axhline(
        half_maximum,
        color="tab:red",
        linestyle=":",
        label="Fitted half-maximum",
    )
    axis.hlines(
        half_maximum,
        left_half_max,
        right_half_max,
        color="tab:red",
        linewidth=2,
    )
    axis.axvline(center, color="tab:purple", linestyle="--", linewidth=1)
    axis.set(
        xlabel="Wavelength (nm; recorded)",
        ylabel="Intensity (a.u.)",
        title="Measured spectrum",
    )
    axis.legend()
    figure.tight_layout()
    figure.savefig(PLOT_DIR / "spectrum.png", dpi=180)
    plt.close(figure)

    return center, fitted_peak_intensity, half_maximum, fwhm, median_intensity


def analyze_autocorrelation(center_wavelength_nm):
    """Fit a tutorial-based IAC model and infer Gaussian pulse duration.

    The model contains the fundamental and second-harmonic fringes described
    by the tutorial's interferometric autocorrelation expression. The dominant
    fundamental fringe has delay period lambda/c and calibrates the recorded
    time axis locally. A Gaussian intensity autocorrelation is sqrt(2) wider
    than the intensity pulse.
    """
    speed_of_light = 299_792_458.0
    optical_fringe_delay = center_wavelength_nm * 1e-9 / speed_of_light
    # vor2 covers a much shorter time window and truncates the autocorrelation
    # wings, so it is shown in the auxiliary plot but excluded from pulse fits.
    accepted_files = ("vor1", "vor3", "nach1", "nach2", "nach3")
    group_pulse_widths_fs = {"vor": [], "nach": []}

    figure, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
    for filename in accepted_files:
        time, signal = np.loadtxt(DATA_DIR / filename, delimiter=",").T
        sample_spacing = np.median(np.diff(time))
        background = np.median(signal)

        # Smooth squared signal to locate the peak; fit the raw trace below.
        smoothed_energy = uniform_filter1d((signal - background) ** 2, size=301)
        initial_center = time[np.argmax(smoothed_energy)]
        local_mask = np.abs(time - initial_center) <= 0.0005
        local_signal = signal[local_mask] - np.mean(signal[local_mask])
        frequencies = np.fft.rfftfreq(local_signal.size, d=sample_spacing)
        power = np.abs(np.fft.rfft(local_signal * np.hanning(local_signal.size))) ** 2
        initial_frequency = frequencies[1 + np.argmax(power[1:])]

        fit_mask = np.abs(time - initial_center) <= 0.002
        initial_guess = (
            background, 0.25, 0.9, 0.0, 0.15, 0.0,
            initial_center, 0.0008, initial_frequency,
        )
        lower_bounds = (
            background - 0.2, -3, -3, -np.pi, -3, -np.pi,
            time.min(), sample_spacing, initial_frequency * 0.8,
        )
        upper_bounds = (
            background + 0.2, 3, 3, np.pi, 3, np.pi,
            time.max(), 0.01, initial_frequency * 1.2,
        )
        parameters, _ = curve_fit(
            interferometric_autocorrelation_model,
            time[fit_mask],
            signal[fit_mask],
            p0=initial_guess,
            bounds=(lower_bounds, upper_bounds),
            maxfev=100_000,
        )
        center_time, sigma_time, carrier_frequency = parameters[6:9]

        autocorrelation_fwhm_time = 2 * np.sqrt(2 * np.log(2)) * abs(sigma_time)
        autocorrelation_fwhm_delay = (
            autocorrelation_fwhm_time * optical_fringe_delay * carrier_frequency
        )
        pulse_width_fs = autocorrelation_fwhm_delay / np.sqrt(2) * 1e15
        group = "vor" if filename.startswith("vor") else "nach"
        group_pulse_widths_fs[group].append(pulse_width_fs)

        relative_delay_fs = (
            (time[fit_mask] - center_time)
            * optical_fringe_delay
            * carrier_frequency
            * 1e15
        )
        axis = axes[0] if group == "vor" else axes[1]
        measured = signal[fit_mask]
        normalization = np.max(np.abs(measured - parameters[0]))
        measured_line, = axis.plot(
            relative_delay_fs,
            (measured - parameters[0]) / normalization,
            linewidth=0.7,
            alpha=0.55,
            label=filename,
        )
        gaussian_envelope = np.exp(
            -0.5 * ((time[fit_mask] - center_time) / sigma_time) ** 2
        )
        envelope_amplitude = (
            parameters[1] + abs(parameters[2]) + abs(parameters[4])
        )
        fitted_envelope = envelope_amplitude * gaussian_envelope
        axis.plot(
            relative_delay_fs,
            np.sqrt(np.clip(fitted_envelope / normalization, 0, None)),
            color=measured_line.get_color(),
            linestyle="--",
            linewidth=2.2,
            label=f"{filename} fit: square-root envelope",
        )
        print(
            f"{filename}: fitted fringe frequency {carrier_frequency / 1e3:.2f} kHz, "
            f"autocorrelation FWHM {autocorrelation_fwhm_delay * 1e15:.1f} fs, "
            f"pulse FWHM {pulse_width_fs:.1f} fs"
        )

    for axis, title in zip(axes, ("Before compressor (vor)", "After compressor (nach)")):
        axis.set_title(title)
        axis.set_xlabel("Relative optical delay (fs; local fringe calibration)")
        axis.grid(True, alpha=0.3)
        axis.legend()
    axes[0].set_ylabel("Normalized signal / square-root fit envelope")
    figure.tight_layout()
    figure.savefig(PLOT_DIR / "autocorrelation_pulsewidth.png", dpi=180)
    plt.close(figure)

    for group, widths in group_pulse_widths_fs.items():
        print(
            f"{group}: pulse duration = {np.mean(widths):.1f} ± "
            f"{np.std(widths, ddof=1):.1f} fs (mean ± sample SD, n={len(widths)})"
        )
    return group_pulse_widths_fs


def main():
    fit_results, slope, intercept, r_squared = analyze_delay_traces()
    (
        peak_wavelength,
        peak_intensity,
        half_maximum,
        spectral_bandwidth,
        median_intensity,
    ) = analyze_spectrum()
    spectral_bandwidth_hz = (
        299_792_458.0 * spectral_bandwidth * 1e-9 / (peak_wavelength * 1e-9) ** 2
    )
    transform_limited_duration_fs = 0.441 / spectral_bandwidth_hz * 1e15

    print(f"Spectrum peak: {peak_wavelength:.2f} nm")
    print(f"Fitted peak intensity: {peak_intensity:.2f} a.u.")
    print(f"Median spectrum intensity: {median_intensity:.2f} a.u.")
    print(f"Fitted half-maximum level: {half_maximum:.2f} a.u.")
    print(f"Gaussian-fit FWHM bandwidth: {spectral_bandwidth:.2f} nm")
    print(f"Approximate frequency bandwidth: {spectral_bandwidth_hz / 1e12:.2f} THz")
    print(f"Gaussian transform-limit duration: {transform_limited_duration_fs:.1f} fs")
    analyze_autocorrelation(peak_wavelength)
    print("\nFile, nominal period, q (m^-1), fitted lifetime (dataset units):")
    for filename, period, wave_number, lifetime in fit_results:
        print(f"{filename}: {period:g}, {wave_number:.6g}, {lifetime:.6g}")

    print("\nExploratory rate-versus-q^2 fit:")
    print(f"  slope = {slope:.6g}")
    print(f"  intercept = {intercept:.6g}")
    print(f"  R^2 = {r_squared:.3f}")
    print("These are not physical transport results without verified calibrations.")


if __name__ == "__main__":
    main()
