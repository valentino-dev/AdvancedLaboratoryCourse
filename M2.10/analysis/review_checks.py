"""Reproduce fit-window, baseline and exclusion sensitivities from the review.

Run after analyze.py. These are sensitivity checks, not additional measurements;
all scan masks are held fixed at the primary 0.5 ps screening threshold.
"""
import json

import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import uniform_filter1d
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
from scipy.stats import linregress

import analyze as analysis


def plot_fringe_calibration(primary):
    """Compare observed peak spacings with the fitted nonuniform scan rate."""
    figure, axes = plt.subplots(2, 1, figsize=(9, 9), constrained_layout=True)
    for axis, (group, name) in zip(axes, primary['selected_autocorrelation'].items()):
        fit = next(fit for fit in primary['autocorrelation'][group] if fit['filename'] == name)
        p = np.array(fit['parameters'])
        time, signal = np.loadtxt(analysis.DATA_DIR / name, delimiter=',').T
        time *= 1e3
        energy = uniform_filter1d((signal - np.median(signal))**2, size=301)
        time -= time[np.argmax(energy)]
        keep = abs(time) <= 2.8
        time, signal = time[keep], signal[keep]
        minimum_distance = max(1, int(0.6 / p[8] / np.median(np.diff(time))))
        peaks, _ = find_peaks(signal, prominence=0.05 * np.ptp(signal), distance=minimum_distance)
        peak_times = time[peaks]
        midpoints = (peak_times[1:] + peak_times[:-1]) / 2
        frequencies = 1 / np.diff(peak_times)
        axis.scatter(midpoints, frequencies, s=14, label=r'Measured: $f_i=1/(x_{i+1}-x_i)$ between maxima')
        relative = time - p[6]
        predicted = p[8] + 2 * p[9] * relative + 3 * p[10] * relative**2
        axis.plot(time, predicted, color='tab:red', linewidth=1.6, label='Red curve: fitted nonlinear scan frequency')
        axis.axhline(p[8], color='0.4', linestyle='--', label='Dashed line: constant-speed approximation')
        axis.set(title=f"{'Before' if group == 'vor' else 'After'} compressor ({name})", xlabel='Recorded time relative to initial peak (ms)',
                 ylabel='Fringe frequency (cycles/ms)')
        axis.text(
            0.5, 0.03,
            'Wider fringe spacing toward the edges means a slower delay scan.\n'
            'The curved trend supports a nonlinear time-to-delay calibration.',
            transform=axis.transAxes, ha='center', va='bottom', fontsize=10,
            bbox=dict(facecolor='white', edgecolor='0.8', alpha=0.95),
        )
        axis.legend(loc='upper center', bbox_to_anchor=(0.5, -0.20), fontsize=10)
    figure.suptitle('Autocorrelation fringe spacing reveals a nonuniform scan', fontsize=14)
    figure.savefig(analysis.PLOT_DIR / 'fringe_calibration.png', dpi=180, bbox_inches='tight')
    plt.close(figure)


def main():
    primary = json.loads((analysis.PROJECT_DIR / 'analysis/results.json').read_text())
    center_nm = primary['spectrum']['center_nm']
    plot_fringe_calibration(primary)
    checks = {'autocorrelation_window': [], 'spectrum_window': [],
              'delay_start': [], 'pump_probe_end': [], 'edge_signal': []}
    for group, name in primary['selected_autocorrelation'].items():
        for window_ms in (2.0, 2.4, 2.8, 3.2):
            fit = analysis.fit_autocorrelation(name, center_nm, window_ms)
            checks['autocorrelation_window'].append({
                'name': name, 'half_window_ms': window_ms,
                'pulse_fs': fit['pulse_width_fs'], 'normalized_rms': fit['normalized_rms'],
            })
        # Same corrected overlap functions, but a constant scan speed.
        fit = analysis.fit_autocorrelation(name, center_nm, nonlinear_scan=False)
        checks.setdefault('constant_scan_speed', []).append({
            'name': name, 'pulse_fs': fit['pulse_width_fs'],
            'normalized_rms': fit['normalized_rms'],
        })
    wavelength, intensity = np.loadtxt(analysis.DATA_DIR / 'spektrum.csv', delimiter=',').T
    for half_width_nm in (None, 50, 30, 20):
        keep = np.ones(wavelength.size, dtype=bool) if half_width_nm is None else abs(wavelength - 789) < half_width_nm
        p, _ = curve_fit(analysis.gaussian_with_offset, wavelength[keep], intensity[keep],
                         p0=(11000, 789, 5, 2500))
        checks['spectrum_window'].append({'half_window_nm': half_width_nm,
                                          'fwhm_nm': analysis.GAUSSIAN_FWHM * abs(p[2])})
    means, q_squared = {}, []
    for name in analysis.TRACE_FILES:
        data = analysis.load_trace(name)
        time, scans = data[:, 0], data[:, 1:]
        _, excluded = analysis.find_anomalous_scans(time, scans)
        means[name] = (time, scans[:, ~excluded].mean(axis=1))
        q_squared.append(analysis.nominal_period_and_wave_number(name)[1]**2)
    original_start = analysis.FIT_START_PS
    try:
        for start_ps in (0.2, 0.5, 1.0, 2.0):
            analysis.FIT_START_PS = start_ps
            lifetimes = []
            for time, mean in means.values():
                p, _ = analysis.fit_trace(time, mean)
                lifetimes.append([p[1], p[1] + p[3]])
            rates = 1 / np.array(lifetimes)[:, 1]
            regression = linregress(q_squared[1:], rates[1:])
            checks['delay_start'].append({'start_ps': start_ps, 'lifetimes_ps': lifetimes,
                                          'diffusivity_m2_s': regression.slope * 1e12})
    finally:
        analysis.FIT_START_PS = original_start
    time, mean = means[analysis.TRACE_FILES[0]]
    for end_ps in (33, 57, 100, 2000):
        keep = time <= end_ps
        p, _ = analysis.fit_trace(time[keep], mean[keep])
        checks['pump_probe_end'].append({'end_ps': end_ps, 'fast_ps': p[1], 'slow_ps': p[1] + p[3]})
    for name in ('vor1', 'vor2', 'vor3'):
        time, signal = np.loadtxt(analysis.DATA_DIR / name, delimiter=',').T
        reference = np.median(signal)
        edge_count = len(signal) // 10
        checks['edge_signal'].append({
            'name': name, 'median': reference,
            'first_10_percent_rms': np.sqrt(np.mean((signal[:edge_count] - reference)**2)),
            'last_10_percent_rms': np.sqrt(np.mean((signal[-edge_count:] - reference)**2)),
        })
    # Show sensitivity to retaining every TG scan, including both flagged repeats.
    rates = []
    for name in analysis.TRACE_FILES[1:]:
        data = analysis.load_trace(name)
        p, _ = analysis.fit_trace(data[:, 0], data[:, 1:].mean(axis=1))
        rates.append(1 / (p[1] + p[3]))
    checks['all_scans_diffusivity_m2_s'] = linregress(q_squared[1:], rates).slope * 1e12
    output = analysis.PROJECT_DIR / 'analysis/review_checks.json'
    output.write_text(json.dumps(checks, indent=2) + '\n')
    print(json.dumps(checks, indent=2))


if __name__ == '__main__':
    main()
