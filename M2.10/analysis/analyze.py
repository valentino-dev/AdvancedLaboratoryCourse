"""Analyze the supplied GaAs transient-grating, spectrum, and IAC data."""

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
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

FIT_START_PS = 0.5
PHASE_MASK_TO_SAMPLE_PERIOD = 0.25
MICROMETRES_TO_METRES = 1e-6
MIN_SCANS_FOR_QUALITY_SCREEN = 8
ROBUST_SCAN_THRESHOLD = 3.5
BOOTSTRAP_REPLICATES = 300

plt.style.use("seaborn-v0_8-whitegrid")
plt.rcParams.update({
    "font.size": 12,
    "axes.labelsize": 12,
    "axes.titlesize": 13,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 10,
})


def load_trace(filename):
    """Load a trace: first column is delay, remaining columns are repeats."""
    return np.loadtxt(DATA_DIR / filename)


def biexponential_with_offset(delay, fast_amplitude, fast_lifetime,
                              slow_amplitude, lifetime_gap, offset):
    """Two decays with ordered lifetimes and a constant signal offset."""
    slow_lifetime = fast_lifetime + lifetime_gap
    return (
        offset
        + fast_amplitude * np.exp(-delay / fast_lifetime)
        + slow_amplitude * np.exp(-delay / slow_lifetime)
    )


def find_anomalous_scans(delay, scans):
    """Flag whole scans whose standardized deviation from the scan median is extreme."""
    score = np.zeros(scans.shape[1])
    excluded = np.zeros(scans.shape[1], dtype=bool)
    if scans.shape[1] < MIN_SCANS_FOR_QUALITY_SCREEN:
        return score, excluded

    fit_region = delay >= FIT_START_PS
    row_median = np.median(scans, axis=1)
    row_mad = np.median(np.abs(scans - row_median[:, None]), axis=1)
    robust_scale = 1.4826 * row_mad[fit_region]
    scale_floor = np.quantile(robust_scale, 0.25)
    robust_scale = np.maximum(robust_scale, scale_floor)
    standardized = (
        scans[fit_region] - row_median[fit_region, None]
    ) / robust_scale[:, None]
    score = np.sqrt(np.mean(standardized**2, axis=0))

    score_median = np.median(score)
    score_mad = np.median(np.abs(score - score_median))
    if score_mad > 0:
        modified_z = 0.6745 * (score - score_median) / score_mad
        excluded = modified_z > ROBUST_SCAN_THRESHOLD
    return score, excluded


def plot_scan_diagnostics(trace_data, output_file):
    """Show retained and excluded complete scans for each delay file."""
    figure, axes = plt.subplots(
        len(trace_data), 1, figsize=(9, 12.5), constrained_layout=False
    )
    for axis, (filename, delay, scans, excluded) in zip(axes, trace_data):
        if scans.shape[1] < MIN_SCANS_FOR_QUALITY_SCREEN:
            label = filename.removeprefix("GaAs_RT_").removesuffix("_X.dat")
            axis.text(
                0.5, 0.5, f"All {scans.shape[1]} scans retained after inspection",
                transform=axis.transAxes, ha="center", va="center",
            )
            axis.set_ylabel(f"{label}\nscan #")
            axis.set_yticks([])
            axis.set_xlim(delay.min(), delay.max())
            axis.grid(True, alpha=0.25)
            continue
        scan_numbers = np.arange(1, scans.shape[1] + 1)
        delay_grid = np.broadcast_to(delay[:, None], scans.shape)
        scan_grid = np.broadcast_to(scan_numbers[None, :], scans.shape)
        axis.scatter(delay_grid[:, ~excluded].ravel(), scan_grid[:, ~excluded].ravel(),
                     s=3, alpha=0.45, color="0.62")
        if excluded.any():
            axis.scatter(delay_grid[:, excluded].ravel(), scan_grid[:, excluded].ravel(),
                         s=14, marker="x", linewidths=1.0, color="tab:red")
        label = filename.removeprefix("GaAs_RT_").removesuffix("_X.dat")
        if "_PP_" not in filename:
            label = f"mask {label} µm"
        axis.set_ylabel(f"{label}\nscan #")
        axis.set_yticks(scan_numbers)
        axis.grid(True, alpha=0.25)
    axes[-1].set_xlabel("Delay (ps)")
    figure.legend(
        handles=(
            Line2D([], [], marker="o", linestyle="", color="0.62", label="Retained scans"),
            Line2D([], [], marker="x", linestyle="", color="tab:red", label="Excluded scan"),
        ),
        loc="upper center",
        bbox_to_anchor=(0.5, 0.965),
        ncol=2,
        fontsize=11,
    )
    figure.suptitle("Whole-scan quality screen", y=0.995)
    figure.subplots_adjust(top=0.92, hspace=0.35, left=0.13, right=0.98)
    figure.savefig(output_file, dpi=180)
    plt.close(figure)


def plot_excluded_scan_comparison(trace_data, output_file):
    """Compare the two flagged scans with the retained population and median."""
    target = next(item for item in trace_data if item[0] == "GaAs_RT_8_X.dat")
    filename, delay, scans, excluded = target
    retained = scans[:, ~excluded]
    median_trace = np.median(scans, axis=1)
    lower, upper = np.quantile(retained, [0.1, 0.9], axis=1)
    scores, _ = find_anomalous_scans(delay, scans)
    score_median = np.median(scores)
    score_mad = np.median(np.abs(scores - score_median))
    modified_scores = 0.6745 * (scores - score_median) / score_mad

    figure, axes = plt.subplots(1, int(excluded.sum()), figsize=(11, 4.4), sharey=True)
    for axis, scan_index in zip(axes, np.flatnonzero(excluded)):
        axis.fill_between(
            delay, lower * 1e6, upper * 1e6, color="0.82", alpha=0.8,
            label="Retained scans: 10th–90th percentile",
        )
        axis.plot(delay, median_trace * 1e6, color="black", linewidth=1.7,
                  label="All-scan pointwise median")
        axis.plot(
            delay, scans[:, scan_index] * 1e6, color="tab:red", linewidth=1.2,
            label=f"Excluded scan {scan_index + 1} (modified score {modified_scores[scan_index]:.2f})",
        )
        axis.axvline(FIT_START_PS, color="0.4", linestyle=":", linewidth=1,
                     label="Fit/screen region starts at 0.5 ps")
        axis.set_title(f"Flagged scan {scan_index + 1}")
        axis.set_xlabel("Delay (ps)")
        axis.grid(True, alpha=0.25)
        axis.legend(fontsize=10, loc="best", frameon=True)
    axes[0].set_ylabel(r"Signal ($10^{-6}$ a.u.)")
    figure.suptitle("Why scans 4 and 10 were flagged (phase-mask period 8 µm)")
    figure.tight_layout()
    figure.savefig(output_file, dpi=180)
    plt.close(figure)


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
    """Fit fast and slow exponential components after time zero."""
    fit_mask = delay >= FIT_START_PS
    fit_delay = delay[fit_mask]
    fit_signal = signal[fit_mask] / 1e-6
    tail = fit_signal[max(0, int(0.8 * fit_signal.size)):]
    offset_guess = float(np.median(tail))
    amplitude_guess = max(float(fit_signal[0] - offset_guess), 1e-3)
    initial_guess = (0.65 * amplitude_guess, 1.0, 0.35 * amplitude_guess, 12.0, offset_guess)
    parameters, covariance = curve_fit(
        biexponential_with_offset,
        fit_delay,
        fit_signal,
        p0=initial_guess,
        bounds=((0.0, 0.05, 0.0, 0.05, -np.inf),
                (np.inf, 500.0, np.inf, 5000.0, np.inf)),
        maxfev=100_000,
    )
    parameters[[0, 2, 4]] *= 1e-6
    return parameters, fit_mask


def nominal_period_and_wave_number(filename):
    """Convert phase-mask period suffix to the actual sample period and wave number."""
    if "_PP_" in filename:
        return 0.0, 0.0

    phase_mask_period_um = float(filename.split("_")[2])
    period_label = phase_mask_period_um * PHASE_MASK_TO_SAMPLE_PERIOD
    period_metres = period_label * MICROMETRES_TO_METRES
    wave_number = 2 * np.pi / period_metres
    return period_label, wave_number


def analyze_delay_traces():
    """Screen scans, fit biexponentials, and regress slow rates against q squared."""
    figure, axes = plt.subplots(1, 2, figsize=(10, 6))
    fit_results = []
    diagnostics = []
    residual_diagnostics = []
    scan_data = {}
    excluded_scan_count = 0
    screened_scan_count = 0

    for filename in TRACE_FILES:
        data = load_trace(filename)
        delay = data[:, 0]
        scans = data[:, 1:]
        scan_data[filename] = (delay, scans)
        quality_scores, excluded = find_anomalous_scans(delay, scans)
        retained_scans = scans[:, ~excluded]
        mean_signal = retained_scans.mean(axis=1)
        standard_error = retained_scans.std(axis=1, ddof=1) / np.sqrt(retained_scans.shape[1])
        diagnostics.append((filename, delay, scans, excluded))
        excluded_scan_count += int(excluded.sum())
        if scans.shape[1] >= MIN_SCANS_FOR_QUALITY_SCREEN:
            screened_scan_count += scans.shape[1]
        excluded_numbers = (np.flatnonzero(excluded) + 1).tolist()
        if scans.shape[1] < MIN_SCANS_FOR_QUALITY_SCREEN:
            print(f"{filename}: inspected {scans.shape[1]} scans; none excluded")
        else:
            print(f"{filename}: excluded whole scan(s) {excluded_numbers or 'none'}")
            print(f"  standardized scan-deviation scores: {np.round(quality_scores, 2)}")
            score_median = np.median(quality_scores)
            score_mad = np.median(np.abs(quality_scores - score_median))
            modified_scores = 0.6745 * (quality_scores - score_median) / score_mad
            print(f"  modified scan scores: {np.round(modified_scores, 2)}")

        parameters, fit_mask = fit_trace(delay, mean_signal)
        fast_amplitude, fast_lifetime, slow_amplitude, lifetime_gap, offset = parameters
        slow_lifetime = fast_lifetime + lifetime_gap
        nominal_period, wave_number = nominal_period_and_wave_number(filename)
        fit_results.append((filename, nominal_period, wave_number, fast_lifetime,
                            slow_lifetime, fast_amplitude, slow_amplitude, offset))
        fit_values = biexponential_with_offset(delay[fit_mask], *parameters)
        fit_residuals = mean_signal[fit_mask] - fit_values
        residual_diagnostics.append((filename, delay[fit_mask], fit_residuals,
                                     np.ptp(mean_signal[fit_mask])))

        label = filename.removeprefix("GaAs_RT_").removesuffix("_X.dat")
        if "_PP_" not in filename:
            label = f"mask {label} µm"
        axis = axes[0] if "_PP_" in filename else axes[1]
        axis.errorbar(
            delay,
            mean_signal,
            yerr=standard_error,
            fmt=".",
            markersize=2.5,
            label=label,
        )
        axis.plot(
            delay[fit_mask],
            biexponential_with_offset(delay[fit_mask], *parameters),
            linewidth=1.5,
        )

    axes[0].set(title="Pump-probe ($q=0$)", xlabel="Delay (ps)")
    axes[1].set(title="Transient grating (finite $q$)", xlabel="Delay (ps)")
    axes[0].set_ylabel("Signal (a.u.)")
    for axis in axes:
        axis.set_title(axis.get_title(), fontsize=15)
        axis.set_xlabel(axis.get_xlabel(), fontsize=14)
        axis.tick_params(axis="both", labelsize=12)
    axes[0].set_ylabel("Signal (a.u.)", fontsize=14)
    figure.suptitle("GaAs pump-probe and transient-grating traces", fontsize=16)
    pp_handles, pp_labels = axes[0].get_legend_handles_labels()
    tg_handles, tg_labels = axes[1].get_legend_handles_labels()
    handles = [*pp_handles, *tg_handles]
    labels = [*pp_labels, *tg_labels]
    handles.append(Line2D([], [], color="0.35", linewidth=1.5))
    labels.append("Solid lines: biexponential fits")
    figure.legend(
        handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.92),
        ncol=4, fontsize=12, frameon=True,
    )
    figure.subplots_adjust(top=0.76, wspace=0.12)
    axes[1].ticklabel_format(axis="y", style="sci", scilimits=(0, 0), useMathText=True)
    figure.savefig(PLOT_DIR / "gaas_decay.png", dpi=180)
    plt.close(figure)

    plot_scan_diagnostics(diagnostics, PLOT_DIR / "outlier_diagnostics.png")
    plot_excluded_scan_comparison(
        diagnostics, PLOT_DIR / "excluded_scan_comparison.png"
    )
    residual_figure, residual_axes = plt.subplots(5, 1, figsize=(8.5, 8), sharex=True)
    for axis, (filename, fit_delay, residuals, signal_range) in zip(
        residual_axes, residual_diagnostics
    ):
        label = filename.removeprefix("GaAs_RT_").removesuffix("_X.dat")
        axis.axhline(0, color="0.25", linewidth=0.8)
        axis.plot(fit_delay, residuals / signal_range, ".", markersize=2.2)
        axis.set_ylabel(f"{label}\n$resid./range$")
        axis.grid(True, alpha=0.25)
    residual_axes[-1].set_xlabel("Delay (ps)")
    residual_figure.suptitle("Biexponential fit residuals (normalized by fitted-window range)")
    residual_figure.tight_layout()
    residual_figure.savefig(PLOT_DIR / "gaas_fit_residuals.png", dpi=180)
    plt.close(residual_figure)
    print(f"Whole scans excluded: {excluded_scan_count}/{screened_scan_count} screened TG scans.")

    wave_numbers = np.array([result[2] for result in fit_results])
    slow_lifetimes = np.array([result[4] for result in fit_results])
    rates = 1 / slow_lifetimes
    q_squared = wave_numbers**2
    q_squared_in_1e12_per_m2 = q_squared * 1e-12
    finite_grating = wave_numbers > 0
    tg_regression = linregress(q_squared[finite_grating], rates[finite_grating])
    tg_slope, tg_intercept, tg_r_value = (
        tg_regression.slope, tg_regression.intercept, tg_regression.rvalue
    )
    all_slope, all_intercept, all_r_value, _, _ = linregress(q_squared, rates)

    # Sensitivity check: retain borderline scan 4 in the mask-period-8 file,
    # while continuing to exclude the clearly divergent scan 10.
    delay_8, scans_8 = scan_data["GaAs_RT_8_X.dat"]
    sensitivity_mask = np.ones(scans_8.shape[1], dtype=bool)
    sensitivity_mask[9] = False
    mean_8_sensitivity = scans_8[:, sensitivity_mask].mean(axis=1)
    parameters_8_sensitivity, _ = fit_trace(delay_8, mean_8_sensitivity)
    slow_lifetime_8_sensitivity = parameters_8_sensitivity[1] + parameters_8_sensitivity[3]
    sensitivity_rates = rates.copy()
    sensitivity_rates[3] = 1 / slow_lifetime_8_sensitivity
    sensitivity_regression = linregress(
        q_squared[finite_grating], sensitivity_rates[finite_grating]
    )

    # Bootstrap whole retained scans to propagate repeat-to-repeat variability
    # through the mean-trace biexponential fits and the finite-q regression.
    rng = np.random.default_rng(20261019)
    bootstrap_diffusivities = []
    bootstrap_lifetimes = []
    for _ in range(BOOTSTRAP_REPLICATES):
        bootstrap_rates = []
        bootstrap_iteration_lifetimes = []
        for filename in TRACE_FILES:
            delay, scans = scan_data[filename]
            _, excluded = find_anomalous_scans(delay, scans)
            retained = scans[:, ~excluded]
            sampled = retained[:, rng.integers(0, retained.shape[1], retained.shape[1])]
            bootstrap_mean = sampled.mean(axis=1)
            bootstrap_parameters, _ = fit_trace(delay, bootstrap_mean)
            fast_lifetime = bootstrap_parameters[1]
            slow_lifetime = bootstrap_parameters[1] + bootstrap_parameters[3]
            bootstrap_iteration_lifetimes.append((fast_lifetime, slow_lifetime))
            bootstrap_rates.append(1 / slow_lifetime)
        bootstrap_lifetimes.append(bootstrap_iteration_lifetimes)
        q_values = np.array([result[2] for result in fit_results])
        bootstrap_fit = linregress(q_values[1:] ** 2, np.asarray(bootstrap_rates)[1:])
        bootstrap_diffusivities.append(bootstrap_fit.slope * 1e12)
    bootstrap_diffusivities = np.asarray(bootstrap_diffusivities)
    figure, axis = plt.subplots(figsize=(7, 4.8))
    axis.scatter(
        q_squared_in_1e12_per_m2[finite_grating], rates[finite_grating],
        label="Finite-$q$ TG slow-component rates",
    )
    axis.scatter(
        q_squared_in_1e12_per_m2[~finite_grating], rates[~finite_grating],
        marker="s", label="$q=0$ pump-probe slow-component rate",
    )
    x_line = np.linspace(0, q_squared_in_1e12_per_m2.max(), 100)
    axis.plot(
        x_line, tg_intercept + tg_slope * x_line * 1e12,
        label="Finite-$q$ TG regression",
    )
    axis.plot(
        x_line, all_intercept + all_slope * x_line * 1e12,
        linestyle="--", label="Regression including $q=0$",
    )
    axis.set(
        xlabel=r"$q^2$ ($10^{12}$ m$^{-2}$)",
        ylabel=r"Slow-component rate (ps$^{-1}$)",
        title="Slow-component rate versus grating wave number squared",
    )
    axis.legend(
        loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=2,
        fontsize=10, frameon=True,
    )
    figure.subplots_adjust(bottom=0.27)
    figure.savefig(PLOT_DIR / "rate_vs_q2.png", dpi=180)
    plt.close(figure)

    return {
        "fits": fit_results,
        "q_squared": q_squared,
        "rates": rates,
        "tg_slope": tg_slope,
        "tg_slope_stderr": tg_regression.stderr,
        "tg_bootstrap_ci": np.percentile(bootstrap_diffusivities, [2.5, 97.5]),
        "lifetime_bootstrap_ci": np.percentile(
            np.asarray(bootstrap_lifetimes), [2.5, 97.5], axis=0
        ),
        "tg_intercept": tg_intercept,
        "tg_r_squared": tg_r_value**2,
        "all_slope": all_slope,
        "all_intercept": all_intercept,
        "all_r_squared": all_r_value**2,
        "scan4_sensitivity_slope": sensitivity_regression.slope,
        "scan4_sensitivity_r_squared": sensitivity_regression.rvalue**2,
    }


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
    parameters, covariance = curve_fit(
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

    figure, axis = plt.subplots(figsize=(8.5, 4.8))
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
    axis.legend(
        loc="center left", bbox_to_anchor=(1.01, 0.5),
        fontsize=10, frameon=True,
    )
    figure.subplots_adjust(right=0.72)
    figure.savefig(PLOT_DIR / "spectrum.png", dpi=180)
    plt.close(figure)

    return center, fitted_peak_intensity, half_maximum, fwhm, median_intensity, covariance


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
    # The truncated vor2 scan does not cover the envelope wings. Fit the
    # remaining candidates, then select the lowest normalized residual per group.
    candidate_files = {"vor": ("vor1", "vor3"), "nach": ("nach1", "nach2", "nach3")}
    candidates = {"vor": [], "nach": []}
    for group, filenames in candidate_files.items():
        for filename in filenames:
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
            parameters, covariance = curve_fit(
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
            pulse_gradient = np.zeros(parameters.size)
            pulse_gradient[7] = pulse_width_fs / sigma_time
            pulse_gradient[8] = pulse_width_fs / carrier_frequency
            pulse_width_stderr_fs = float(
                np.sqrt(pulse_gradient @ covariance @ pulse_gradient)
            )
            fitted = interferometric_autocorrelation_model(time[fit_mask], *parameters)
            residual_fraction = np.sqrt(np.mean((signal[fit_mask] - fitted) ** 2)) / np.ptp(signal[fit_mask])
            candidates[group].append((residual_fraction, filename, time, signal, fit_mask,
                                      parameters, pulse_width_fs, pulse_width_stderr_fs,
                                      autocorrelation_fwhm_delay))
  
    selected = {group: min(group_candidates, key=lambda candidate: candidate[0])
                for group, group_candidates in candidates.items()}
    figure, axes = plt.subplots(2, 1, figsize=(9, 9.5), sharex=True)
    selected_results = {}
    for axis, (group, title) in zip(axes, (("vor", "Before compressor"), ("nach", "After compressor"))):
        (residual_fraction, filename, time, signal, fit_mask, parameters,
         pulse_width_fs, pulse_width_stderr_fs, autocorrelation_fwhm_delay) = selected[group]
        center_time, sigma_time, carrier_frequency = parameters[6:9]
        selected_results[group] = (
            filename, pulse_width_fs, pulse_width_stderr_fs, residual_fraction
        )
        relative_delay_fs = ((time[fit_mask] - center_time) * optical_fringe_delay
                             * carrier_frequency * 1e15)
        measured = signal[fit_mask]
        normalization = np.max(np.abs(measured - parameters[0]))
        measured_line, = axis.plot(
            relative_delay_fs, (measured - parameters[0]) / normalization,
            linewidth=0.55, alpha=0.45, label=f"Selected trace: {filename}",
        )
        gaussian_envelope = np.exp(-0.5 * ((time[fit_mask] - center_time) / sigma_time) ** 2)
        envelope_amplitude = parameters[1] + abs(parameters[2]) + abs(parameters[4])
        fitted_envelope = envelope_amplitude * gaussian_envelope
        axis.plot(
            relative_delay_fs,
            np.sqrt(np.clip(fitted_envelope / normalization, 0, None)),
            color=measured_line.get_color(), linestyle="--", linewidth=1.8,
            label="Fit: square-root envelope",
        )
        axis.set_title(title, fontsize=15)
        axis.set_xlabel("Relative optical delay (fs; local fringe calibration)", fontsize=14)
        axis.tick_params(axis="both", labelsize=12)
        axis.grid(True, alpha=0.3)
        axis.legend(loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=1,
                    fontsize=11, frameon=True)
        print(
            f"{group}: selected {filename} (normalized fit residual {residual_fraction:.4f}); "
            f"fringe frequency {carrier_frequency / 1e3:.2f} kHz, "
            f"pulse FWHM {pulse_width_fs:.1f} fs"
        )
    for group, group_candidates in candidates.items():
        for candidate in group_candidates:
            residual_fraction, filename, *_, pulse_width_fs, pulse_width_stderr_fs, _ = candidate
            print(
                f"  candidate {filename}: {pulse_width_fs:.2f} +/- "
                f"{pulse_width_stderr_fs:.2f} fs (fit covariance), "
                f"normalized residual {residual_fraction:.4f}"
            )

    for axis in axes:
        axis.set_ylabel("Normalized signal and fitted envelope", fontsize=14)
    figure.subplots_adjust(left=0.12, right=0.995, bottom=0.08, top=0.98, hspace=0.62)
    figure.savefig(PLOT_DIR / "autocorrelation_pulsewidth.png", dpi=180)
    plt.close(figure)

    return selected_results


def main():
    delay_results = analyze_delay_traces()
    (
        peak_wavelength,
        peak_intensity,
        half_maximum,
        spectral_bandwidth,
        median_intensity,
        spectrum_covariance,
    ) = analyze_spectrum()
    spectral_bandwidth_hz = (
        299_792_458.0 * spectral_bandwidth * 1e-9 / (peak_wavelength * 1e-9) ** 2
    )
    transform_limited_duration_fs = 0.441 / spectral_bandwidth_hz * 1e15

    print(f"Spectrum peak: {peak_wavelength:.2f} nm")
    print(f"Fitted peak intensity: {peak_intensity:.2f} a.u.")
    print(f"Median spectrum intensity: {median_intensity:.2f} a.u.")
    print(f"Fitted half-maximum level: {half_maximum:.2f} a.u.")
    sigma_error = np.sqrt(spectrum_covariance[2, 2])
    bandwidth_error = 2 * np.sqrt(2 * np.log(2)) * sigma_error
    print(f"Gaussian-fit FWHM bandwidth: {spectral_bandwidth:.2f} +/- {bandwidth_error:.2f} nm (fit covariance)")
    print(f"Approximate frequency bandwidth: {spectral_bandwidth_hz / 1e12:.2f} THz")
    print(f"Gaussian transform-limit duration: {transform_limited_duration_fs:.1f} fs")
    selected_autocorrelations = analyze_autocorrelation(peak_wavelength)
    print("\nFile, sample period (um), q (m^-1), fast lifetime (ps), slow lifetime (ps):")
    lifetime_ci = delay_results["lifetime_bootstrap_ci"]
    for index, (filename, period, wave_number, fast_lifetime, slow_lifetime, *_) in enumerate(delay_results["fits"]):
        fast_ci = lifetime_ci[:, index, 0]
        slow_ci = lifetime_ci[:, index, 1]
        print(
            f"{filename}: {period:g}, {wave_number:.6g}, "
            f"tau_fast={fast_lifetime:.4g} [{fast_ci[0]:.4g}, {fast_ci[1]:.4g}], "
            f"tau_slow={slow_lifetime:.4g} [{slow_ci[0]:.4g}, {slow_ci[1]:.4g}] ps"
        )
    print("\nSelected autocorrelation scans (single cleanest trace per group):")
    for group, (filename, pulse_width_fs, pulse_error_fs, residual) in selected_autocorrelations.items():
        print(
            f"{group}: {filename}, {pulse_width_fs:.2f} +/- {pulse_error_fs:.2f} fs "
            f"(fit covariance), normalized residual {residual:.4f}"
        )

    tg_slope = delay_results["tg_slope"]
    diffusion_m2_s = tg_slope * 1e12
    diffusion_stderr = delay_results["tg_slope_stderr"] * 1e12
    print("\nSlow-component rate versus q^2 regression:")
    print(f"  finite-q TG slope = {tg_slope:.6g} m^2/ps")
    print(f"  corresponding D = ({diffusion_m2_s:.4g} +/- {diffusion_stderr:.2g}) m^2/s")
    print(
        "  repeat-scan bootstrap 95% interval = "
        f"[{delay_results['tg_bootstrap_ci'][0]:.4g}, "
        f"{delay_results['tg_bootstrap_ci'][1]:.4g}] m^2/s"
    )
    print(f"  intercept = {delay_results['tg_intercept']:.6g} ps^-1")
    print(f"  TG-only R^2 = {delay_results['tg_r_squared']:.3f}")
    print(
        "  retaining borderline scan 4 (still rejecting scan 10): "
        f"D = {delay_results['scan4_sensitivity_slope'] * 1e12:.4g} m^2/s, "
        f"R^2 = {delay_results['scan4_sensitivity_r_squared']:.3f}"
    )
    print("  including q=0:")
    print(f"    slope = {delay_results['all_slope']:.6g} m^2/ps")
    print(f"    intercept = {delay_results['all_intercept']:.6g} ps^-1")
    print(f"    R^2 = {delay_results['all_r_squared']:.3f}")


if __name__ == "__main__":
    main()
