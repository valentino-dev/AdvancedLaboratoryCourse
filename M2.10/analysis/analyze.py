"""Analyze the supplied GaAs transient-grating, spectrum, and IAC data."""

import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import linregress
from scipy.ndimage import uniform_filter1d
from scipy.signal import find_peaks


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
SPEED_OF_LIGHT = 299_792_458.0
GAUSSIAN_FWHM = 2 * np.sqrt(2 * np.log(2))
GAUSSIAN_TIME_BANDWIDTH = 2 * np.log(2) / np.pi

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


def modified_scan_scores(scores):
    """Median/MAD scores; a zero MAD gives no defensible automatic exclusions."""
    median = np.median(scores)
    mad = np.median(np.abs(scores - median))
    if mad == 0:
        return np.zeros_like(scores)
    return 0.6745 * (scores - median) / mad


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
    scale_floor = max(
        np.quantile(robust_scale, 0.25),
        np.finfo(float).eps * max(float(np.max(np.abs(scans))), np.finfo(float).tiny),
    )
    robust_scale = np.maximum(robust_scale, scale_floor)
    standardized = (
        scans[fit_region] - row_median[fit_region, None]
    ) / robust_scale[:, None]
    score = np.sqrt(np.mean(standardized**2, axis=0))

    excluded = modified_scan_scores(score) > ROBUST_SCAN_THRESHOLD
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
                0.5, 0.5, f"All {scans.shape[1]} scans retained; too few for robust screen",
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
    modified_scores = modified_scan_scores(scores)

    if not excluded.any():
        return
    figure, axes_grid = plt.subplots(
        1, int(excluded.sum()), figsize=(11, 4.4), sharey=True, squeeze=False
    )
    axes = axes_grid.ravel()
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


def fringe_cycles(time_ms, center_ms, frequency_per_ms, quadratic, cubic):
    """Optical-delay coordinate in carrier cycles, allowing a nonlinear mirror scan."""
    relative_time = time_ms - center_ms
    return (frequency_per_ms * relative_time + quadratic * relative_time**2
            + cubic * relative_time**3)


def iac_basis(time_ms, center_ms, pulse_sigma_cycles, frequency_per_ms,
              quadratic, cubic):
    """Gaussian overlap integrals: DC/2ω and ω have different Gaussian widths."""
    cycles = fringe_cycles(time_ms, center_ms, frequency_per_ms, quadratic, cubic)
    intensity_overlap = np.exp(-cycles**2 / (4 * pulse_sigma_cycles**2))
    field_overlap = np.exp(-3 * cycles**2 / (16 * pulse_sigma_cycles**2))
    phase = 2 * np.pi * cycles
    return np.column_stack((
        np.ones_like(cycles), intensity_overlap,
        field_overlap * np.cos(phase), field_overlap * np.sin(phase),
        intensity_overlap * np.cos(2 * phase), intensity_overlap * np.sin(2 * phase),
    ))


def interferometric_autocorrelation_model(time_ms, offset, dc_amplitude,
                                         cos_omega, sin_omega, cos_twice_omega,
                                         sin_twice_omega, center_ms,
                                         pulse_sigma_cycles, frequency_per_ms,
                                         quadratic, cubic):
    """Fit the tutorial overlap functions using sine/cosine amplitudes.

    The fitted pulse sigma is in optical cycles, not the width of an IAC
    envelope. Sine/cosine coefficients avoid artificial phase-boundary minima.
    The cubic recorded-time mapping accounts for changing mirror velocity.
    """
    amplitudes = (offset, dc_amplitude, cos_omega, sin_omega,
                  cos_twice_omega, sin_twice_omega)
    return iac_basis(time_ms, center_ms, pulse_sigma_cycles, frequency_per_ms,
                     quadratic, cubic) @ amplitudes


def fit_trace(delay, signal, return_covariance=False):
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
    scale = np.array([1e-6, 1, 1e-6, 1, 1e-6])
    parameters *= scale
    covariance *= np.outer(scale, scale)
    if return_covariance:
        return parameters, fit_mask, covariance
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
            print(f"{filename}: {scans.shape[1]} scans; robust screen not applied; all retained")
        else:
            print(f"{filename}: excluded whole scan(s) {excluded_numbers or 'none'}")
            print(f"  standardized scan-deviation scores: {np.round(quality_scores, 2)}")
            modified_scores = modified_scan_scores(quality_scores)
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
        points = axis.errorbar(
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
            color=points.lines[0].get_color(),
        )

    axes[0].set(title="Pump-probe ($q=0$; early-time view)", xlabel="Delay (ps)")
    axes[0].set_xlim(-1, 35)
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
    residual_figure, residual_axes = plt.subplots(5, 1, figsize=(9, 11), sharex=False)
    for axis, (filename, fit_delay, residuals, signal_range) in zip(
        residual_axes, residual_diagnostics
    ):
        label = filename.removeprefix("GaAs_RT_").removesuffix("_X.dat")
        axis.axhline(0, color="0.25", linewidth=0.8)
        axis.plot(fit_delay, residuals / signal_range, ".", markersize=2.2)
        axis.set_ylabel(f"{label}\n$resid./range$")
        axis.set_xlabel("Delay (ps)")
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
    bootstrap_rate_ci = np.percentile(
        1 / np.asarray(bootstrap_lifetimes)[:, :, 1], [2.5, 97.5], axis=0
    )
    rate_errors = np.array([rates - bootstrap_rate_ci[0], bootstrap_rate_ci[1] - rates])
    figure, axis = plt.subplots(figsize=(7, 4.8))
    axis.errorbar(
        q_squared_in_1e12_per_m2[finite_grating], rates[finite_grating],
        yerr=rate_errors[:, finite_grating], fmt="o", capsize=3,
        label="Finite-$q$ TG slow-component rates",
    )
    axis.errorbar(
        q_squared_in_1e12_per_m2[~finite_grating], rates[~finite_grating],
        yerr=rate_errors[:, ~finite_grating], fmt="s", capsize=3, label="$q=0$ pump-probe slow-component rate",
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
    axis.plot(
        x_line, sensitivity_regression.intercept + sensitivity_regression.slope * x_line * 1e12,
        linestyle=":", color="tab:purple", label="TG fit retaining borderline scan 4",
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
    figure.subplots_adjust(bottom=0.32)
    figure.savefig(PLOT_DIR / "rate_vs_q2.png", dpi=180, bbox_inches="tight")
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

    figure, axes = plt.subplots(2, 1, figsize=(9, 7.5), constrained_layout=True)
    for axis in axes:
        axis.plot(wavelength, intensity, linewidth=1.2, label="Measured spectrum")
        axis.plot(wavelength, gaussian_with_offset(wavelength, *parameters),
                  linewidth=1.8, label="Gaussian fit")
        axis.axhline(median_intensity, color="tab:green", linestyle="--",
                     label="Whole-record median (descriptive)")
        axis.axhline(half_maximum, color="tab:red", linestyle=":", label="Fitted half-maximum")
        axis.hlines(half_maximum, left_half_max, right_half_max, color="tab:red", linewidth=3)
        axis.set(xlabel="Wavelength (nm)", ylabel="Intensity (a.u.)")
        axis.tick_params(labelsize=11)
    axes[0].set_title("Recorded spectrum: full range")
    axes[1].set(xlim=(center - 25, center + 25), title=f"FWHM = {fwhm:.3f} ± {GAUSSIAN_FWHM * np.sqrt(covariance[2, 2]):.3f} nm (fit SE)")
    axes[1].legend(loc="upper center", bbox_to_anchor=(0.5, -0.21), ncol=2, fontsize=11)
    figure.savefig(PLOT_DIR / "spectrum.png", dpi=180, bbox_inches="tight")
    plt.close(figure)

    return center, fitted_peak_intensity, half_maximum, fwhm, median_intensity, covariance


def fit_autocorrelation(filename, center_wavelength_nm, window_ms=2.8,
                        nonlinear_scan=True):
    """Fit raw IAC samples; initialize scan calibration from successive maxima."""
    time_seconds, signal = np.loadtxt(DATA_DIR / filename, delimiter=",").T
    time_ms = time_seconds * 1e3
    sample_spacing_ms = float(np.median(np.diff(time_ms)))
    background = np.median(signal)
    smoothed_energy = uniform_filter1d((signal - background)**2, size=301)
    initial_center_ms = time_ms[np.argmax(smoothed_energy)]
    fit_mask = np.abs(time_ms - initial_center_ms) <= window_ms
    # Center the numerical coordinate to keep the polynomial well conditioned.
    fit_time = time_ms[fit_mask] - initial_center_ms
    fit_signal = signal[fit_mask]
    local_mask = np.abs(fit_time) <= 0.5
    local_signal = fit_signal[local_mask] - np.mean(fit_signal[local_mask])
    frequencies = np.fft.rfftfreq(local_signal.size, d=sample_spacing_ms)
    power = np.abs(np.fft.rfft(local_signal * np.hanning(local_signal.size)))**2
    initial_frequency = frequencies[1 + np.argmax(power[1:])]
    # Maxima initialize the fringe phase only; every raw sample enters the fit.
    stride = 5
    coarse_time, coarse_signal = fit_time[::stride], fit_signal[::stride]
    minimum_distance = max(1, int(0.6 / initial_frequency / sample_spacing_ms / stride))
    peaks, _ = find_peaks(coarse_signal, prominence=0.05 * np.ptp(coarse_signal),
                          distance=minimum_distance)
    if peaks.size < 10:
        raise ValueError(f"Too few resolved fringes in {filename}")
    cubic, quadratic, frequency, _ = np.polyfit(coarse_time[peaks],
                                               np.arange(peaks.size), 3)
    if not nonlinear_scan:
        quadratic = cubic = 0.0
    shape_guess = (0.0, 0.55 * frequency, frequency, quadratic, cubic)
    basis = iac_basis(fit_time, *shape_guess)
    amplitude_guess = np.linalg.lstsq(basis, fit_signal, rcond=None)[0]
    initial_guess = np.r_[amplitude_guess, shape_guess]
    lower = [-np.inf] * 6 + [-0.5, 0.1, initial_frequency * 0.8, -10, -5]
    upper = [np.inf] * 6 + [0.5, 200, initial_frequency * 1.2, 10, 5]
    if nonlinear_scan:
        fit_model = interferometric_autocorrelation_model
    else:
        def fit_model(x, *parameters):
            return interferometric_autocorrelation_model(x, *parameters, 0.0, 0.0)
        initial_guess = initial_guess[:-2]
        lower, upper = lower[:-2], upper[:-2]
    parameters, covariance = curve_fit(
        fit_model, fit_time, fit_signal, p0=initial_guess,
        bounds=(lower, upper), x_scale="jac", maxfev=20_000,
    )
    if not nonlinear_scan:
        parameters = np.r_[parameters, 0.0, 0.0]
    relative_time = fit_time - parameters[6]
    scan_speed = (parameters[8] + 2 * parameters[9] * relative_time
                  + 3 * parameters[10] * relative_time**2)
    if np.any(scan_speed <= 0):
        raise ValueError(f"Non-monotonic fitted delay mapping in {filename}")
    delay_per_cycle_fs = center_wavelength_nm * 1e-9 / SPEED_OF_LIGHT * 1e15
    pulse_width_fs = GAUSSIAN_FWHM * parameters[7] * delay_per_cycle_fs
    pulse_stderr_fs = GAUSSIAN_FWHM * np.sqrt(covariance[7, 7]) * delay_per_cycle_fs
    fitted_signal = interferometric_autocorrelation_model(fit_time, *parameters)
    normalized_rms = np.sqrt(np.mean((fit_signal - fitted_signal)**2)) / np.ptp(fit_signal)
    return {
        "filename": filename, "parameters": parameters,
        "pulse_width_fs": pulse_width_fs, "pulse_stderr_fs": pulse_stderr_fs,
        "normalized_rms": normalized_rms, "time_ms": fit_time,
        "signal": fit_signal, "fitted_signal": fitted_signal,
        "delay_per_cycle_fs": delay_per_cycle_fs,
        "scan_frequency_range_per_ms": [float(scan_speed.min()), float(scan_speed.max())],
    }


def analyze_autocorrelation(center_wavelength_nm):
    """Fit candidates and display the normalized IAC modulation envelope."""
    # vor2 remains excluded: a separate edge-signal check demonstrates that its
    # record ends before the wings reach baseline; sampling interval alone is not a reason.
    filenames = {"vor": ("vor1", "vor3"), "nach": ("nach1", "nach2", "nach3")}
    candidates = {group: [fit_autocorrelation(name, center_wavelength_nm)
                          for name in names] for group, names in filenames.items()}
    selected = {group: min(fits, key=lambda fit: fit["normalized_rms"])
                for group, fits in candidates.items()}
    figure, axes = plt.subplots(2, 1, figsize=(9, 10), constrained_layout=True)
    for axis, group, title in zip(axes, ("vor", "nach"),
                                  ("Before compressor", "After compressor")):
        fit = selected[group]
        p = fit["parameters"]
        cycles = fringe_cycles(fit["time_ms"], p[6], p[8], p[9], p[10])
        optical_delay_fs = cycles * fit["delay_per_cycle_fs"]
        normalization = np.max(np.abs(fit["signal"] - p[0]))
        axis.plot(optical_delay_fs, (fit["signal"] - p[0]) / normalization,
                  linewidth=0.55, alpha=0.55, label=f"Raw signal: {fit['filename']}")
        g = np.exp(-cycles**2 / (4 * p[7]**2))
        h = np.exp(-3 * cycles**2 / (16 * p[7]**2))
        bound = (abs(p[1]) + np.hypot(p[4], p[5])) * g + np.hypot(p[2], p[3]) * h
        axis.plot(optical_delay_fs, bound / normalization,
                  color="tab:red", linestyle="--", linewidth=2,
                  label=r"Fitted envelope: $U/N$")
        axis.set_title(
            f"{title}: {fit['pulse_width_fs']:.2f} ± {fit['pulse_stderr_fs']:.2f} fs (fit SE)",
            fontsize=14,
        )
        axis.set_xlabel("Optical delay (fs; nonlinear fringe calibration)", fontsize=14)
        axis.set_ylabel("Dimensionless normalized signal", fontsize=14)
        axis.tick_params(axis="both", labelsize=12)
        axis.legend(loc="upper center", bbox_to_anchor=(0.5, -0.18), fontsize=11)
    figure.savefig(PLOT_DIR / "autocorrelation_pulsewidth.png", dpi=180, bbox_inches="tight")
    plt.close(figure)

    # Residuals show model quality without drawing fitted oscillations over the data.
    figure, axes = plt.subplots(2, 1, figsize=(9, 6), constrained_layout=True)
    for axis, group in zip(axes, ("vor", "nach")):
        fit = selected[group]
        residual = (fit["signal"] - fit["fitted_signal"]) / np.ptp(fit["signal"])
        axis.plot(fit["time_ms"], residual, linewidth=0.45)
        axis.axhline(0, color="black", linewidth=0.7)
        axis.set(title=f"{fit['filename']}: normalized RMS {fit['normalized_rms']:.3f}",
                 xlabel="Recorded time relative to initial peak (ms)", ylabel="Residual / range")
    figure.savefig(PLOT_DIR / "autocorrelation_residuals.png", dpi=180)
    plt.close(figure)

    for fits in candidates.values():
        for fit in fits:
            print(f"IAC {fit['filename']}: {fit['pulse_width_fs']:.3f} +/- "
                  f"{fit['pulse_stderr_fs']:.3f} fs (covariance), "
                  f"normalized RMS {fit['normalized_rms']:.5f}, "
                  f"scan frequency range {fit['scan_frequency_range_per_ms']} cycles/ms")
    compact = {group: [{key: value.tolist() if isinstance(value, np.ndarray) else value
                        for key, value in fit.items()
                        if key not in ("time_ms", "signal", "fitted_signal")}
                       for fit in fits] for group, fits in candidates.items()}
    rows = [
        r"% Generated by analysis/analyze.py; do not edit numerical values here.",
        r"\begin{table}[htbp]\centering",
        r"\caption{Gaussian pulse-intensity FWHM from the raw IAC fits. Errors are conditional covariance standard errors at fixed wavelength, excluding calibration and model uncertainty. Normalized RMS is residual RMS divided by the measured range; its uncertainty has not been estimated.}\label{tab:iac-fits}",
        r"\begin{tabular}{lrr}\toprule",
        r"Trace & Pulse FWHM (fs) & Normalized RMS\\\midrule",
    ]
    for group, fits in candidates.items():
        for fit in fits:
            selected_marker = r" (selected)" if fit is selected[group] else ""
            rows.append(
                f"{fit['filename']}{selected_marker} & "
                f"${fit['pulse_width_fs']:.2f} \\pm {fit['pulse_stderr_fs']:.2f}$ & "
                f"{fit['normalized_rms']:.4f}" + r"\\"
            )
    rows.extend([r"\bottomrule\end{tabular}", r"\end{table}"])
    (PROJECT_DIR / "latex" / "autocorrelation_results.tex").write_text("\n".join(rows) + "\n")
    return selected, compact


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
    transform_limited_duration_fs = GAUSSIAN_TIME_BANDWIDTH / spectral_bandwidth_hz * 1e15

    print(f"Spectrum peak: {peak_wavelength:.2f} nm")
    print(f"Fitted peak intensity: {peak_intensity:.2f} a.u.")
    print(f"Median spectrum intensity: {median_intensity:.2f} a.u.")
    print(f"Fitted half-maximum level: {half_maximum:.2f} a.u.")
    sigma_error = np.sqrt(spectrum_covariance[2, 2])
    bandwidth_error = 2 * np.sqrt(2 * np.log(2)) * sigma_error
    print(f"Gaussian-fit FWHM bandwidth: {spectral_bandwidth:.2f} +/- {bandwidth_error:.2f} nm (fit covariance)")
    print(f"Approximate frequency bandwidth: {spectral_bandwidth_hz / 1e12:.2f} THz")
    print(f"Gaussian transform-limit duration: {transform_limited_duration_fs:.1f} fs")
    selected_autocorrelations, autocorrelation_candidates = analyze_autocorrelation(peak_wavelength)
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
    print("\nSelected autocorrelation scans (lowest normalized RMS per group):")
    for group, fit in selected_autocorrelations.items():
        filename = fit["filename"]
        pulse_width_fs, pulse_error_fs = fit["pulse_width_fs"], fit["pulse_stderr_fs"]
        residual = fit["normalized_rms"]
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

    results = {
        "spectrum": {"center_nm": peak_wavelength, "fwhm_nm": spectral_bandwidth,
                     "fwhm_covariance_stderr_nm": bandwidth_error,
                     "frequency_bandwidth_hz": spectral_bandwidth_hz,
                     "transform_limit_fs": transform_limited_duration_fs},
        "delay": delay_results, "autocorrelation": autocorrelation_candidates,
        "selected_autocorrelation": {group: fit["filename"]
                                      for group, fit in selected_autocorrelations.items()},
    }
    def numpy_to_json(value):
        if isinstance(value, np.ndarray):
            return value.tolist()
        if isinstance(value, np.generic):
            return value.item()
        raise TypeError(f"Cannot serialize {type(value)}")
    (PROJECT_DIR / "analysis" / "results.json").write_text(
        json.dumps(results, indent=2, default=numpy_to_json) + "\n"
    )


if __name__ == "__main__":
    main()
