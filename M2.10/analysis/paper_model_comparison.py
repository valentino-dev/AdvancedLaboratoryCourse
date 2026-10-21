"""Evaluate the shared-diffusion assumption in Chen et al., arXiv:1805.01561.

The paper's total heterodyne signal includes a population contribution. The
course data are treated as phase-separated grating amplitudes, so only the
coefficient of cos(heterodyne phase) is compared with the grating rate law.
The extra PP exponential is an empirical diagnostic, not a process assignment.
"""
import json

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import linregress

import analyze as analysis


REPLICATES = 300


def three_exponentials(delay, amplitude_fast, lifetime_fast, amplitude_middle,
                       gap_middle, amplitude_long, gap_long, offset):
    lifetime_middle = lifetime_fast + gap_middle
    lifetime_long = lifetime_middle + gap_long
    return (offset + amplitude_fast * np.exp(-delay / lifetime_fast)
            + amplitude_middle * np.exp(-delay / lifetime_middle)
            + amplitude_long * np.exp(-delay / lifetime_long))


def fit_three_exponentials(delay, signal):
    fit_mask = delay >= analysis.FIT_START_PS
    parameters, covariance = curve_fit(
        three_exponentials, delay[fit_mask], signal[fit_mask] * 1e6,
        p0=(180, 0.5, 40, 3, 10, 400, 30),
        bounds=((0, .05, 0, .05, 0, .05, -np.inf),
                (np.inf, 500, np.inf, 5000, np.inf, 1e5, np.inf)),
        maxfev=100000,
    )
    lifetime_gradient = np.zeros((3, 7))
    lifetime_gradient[0, 1] = 1
    lifetime_gradient[1, [1, 3]] = 1
    lifetime_gradient[2, [1, 3, 5]] = 1
    lifetimes = lifetime_gradient @ parameters
    errors = np.sqrt(np.diag(lifetime_gradient @ covariance @ lifetime_gradient.T))
    return parameters, lifetimes, errors, fit_mask


def aicc(residuals, parameter_count):
    """Gaussian unknown common variance; include that variance in parameter count."""
    count = residuals.size
    k = parameter_count + 1
    return (count * np.log(np.mean(residuals**2)) + 2*k
            + 2*k*(k+1)/(count-k-1))


def main():
    data = analysis.load_trace('GaAs_RT_PP_X.dat')
    delay, scans = data[:, 0], data[:, 1:]
    mean_signal = scans.mean(axis=1)
    parameters_two, mask = analysis.fit_trace(delay, mean_signal)
    parameters_three, lifetimes, errors, _ = fit_three_exponentials(delay, mean_signal)
    two_fit = analysis.biexponential_with_offset(delay[mask], *parameters_two)
    three_fit = three_exponentials(delay[mask], *parameters_three) * 1e-6
    normalization = np.ptp(mean_signal[mask])
    two_residuals = (mean_signal[mask] - two_fit) / normalization
    three_residuals = (mean_signal[mask] - three_fit) / normalization
    results = {
        'paper': 'https://arxiv.org/abs/1805.01561',
        'pp_three_exponential_lifetimes_ps': lifetimes.tolist(),
        'pp_three_exponential_covariance_SE_ps': errors.tolist(),
        'pp_normalized_RMS_two': float(np.sqrt(np.mean(two_residuals**2))),
        'pp_normalized_RMS_three': float(np.sqrt(np.mean(three_residuals**2))),
        'pp_AICc_two': float(aicc(two_residuals, 5)),
        'pp_AICc_three': float(aicc(three_residuals, 7)),
    }
    repeat_lifetimes = [fit_three_exponentials(delay, scans[:, index])[1].tolist()
                        for index in range(scans.shape[1])]
    results['pp_three_exponential_repeat_lifetimes_ps'] = repeat_lifetimes

    q_squared = []
    rates = []
    retained = []
    for filename in analysis.TRACE_FILES[1:]:
        trace = analysis.load_trace(filename)
        times, repeats = trace[:, 0], trace[:, 1:]
        _, excluded = analysis.find_anomalous_scans(times, repeats)
        repeats = repeats[:, ~excluded]
        retained.append((times, repeats))
        parameters, _ = analysis.fit_trace(times, repeats.mean(axis=1))
        rates.append((1/parameters[1], 1/(parameters[1]+parameters[3])))
        q_squared.append(analysis.nominal_period_and_wave_number(filename)[1]**2)
    q_squared, rates = np.array(q_squared), np.array(rates)
    regressions = [linregress(q_squared, rates[:, component]) for component in (0, 1)]
    for name, regression in zip(('fast', 'slow'), regressions):
        results[name+'_diffusion_m2_s'] = regression.slope * 1e12
        results[name+'_diffusion_regression_SE_m2_s'] = regression.stderr * 1e12

    # Paired bootstrap preserves the correlation between fast/slow components
    # obtained from fitting the same resampled complete scans.
    rng = np.random.default_rng(20261020)
    bootstrap_slopes = []
    bootstrap_rates = []
    for _ in range(REPLICATES):
        iteration_rates = []
        for times, repeats in retained:
            columns = rng.integers(0, repeats.shape[1], repeats.shape[1])
            parameters, _ = analysis.fit_trace(times, repeats[:, columns].mean(axis=1))
            iteration_rates.append((1/parameters[1], 1/(parameters[1]+parameters[3])))
        iteration_rates = np.array(iteration_rates)
        bootstrap_rates.append(iteration_rates)
        bootstrap_slopes.append([linregress(q_squared, iteration_rates[:, component]).slope*1e12
                                 for component in (0, 1)])
    bootstrap_slopes = np.array(bootstrap_slopes)
    rate_intervals = np.percentile(np.array(bootstrap_rates), [2.5, 97.5], axis=0)
    results['component_rate_CI95_per_ps'] = rate_intervals.tolist()
    results['paired_bootstrap_fast_minus_slow_diffusion_CI95_m2_s'] = np.percentile(
        bootstrap_slopes[:, 0]-bootstrap_slopes[:, 1], [2.5, 97.5]).tolist()
    results['paired_bootstrap_component_diffusion_CI95_m2_s'] = np.percentile(
        bootstrap_slopes, [2.5, 97.5], axis=0).tolist()
    results['uncertainty_scope'] = (
        'SE denotes conditional covariance or regression standard error. Bootstrap '
        'CI resamples retained scans with fixed screening. Supplied delays and grating '
        'spacings are treated as exact. Acquisition drift, model and selection uncertainties are unquantified. RMS and AICc '
        'are descriptive diagnostics without estimated uncertainties. AICc assumes '
        'independent equal-variance Gaussian errors; PP residual structure limits it.'
    )

    figure, axes = plt.subplots(2, 1, figsize=(9, 8), constrained_layout=True)
    axes[0].plot(delay[mask], two_residuals, '.-', markersize=3, linewidth=.8,
                 label='Two exponentials + constant')
    axes[0].plot(delay[mask], three_residuals, '.-', markersize=3, linewidth=.8,
                 label='Three exponentials + constant (diagnostic)')
    axes[0].axhline(0, color='0.3', linewidth=.7)
    axes[0].set(xlabel='Delay (ps)', ylabel='Residual / measured range',
                title='Pump–probe: testing an additional slow contribution')
    axes[0].legend(fontsize=11)
    grid = np.linspace(0, q_squared.max(), 100)
    for component, (name, regression) in enumerate(zip(('Fast', 'Slow'), regressions)):
        line, = axes[1].plot(grid*1e-12, regression.intercept + regression.slope*grid,
                            label=f'{name} rate regression')
        errors = np.array([rates[:, component]-rate_intervals[0, :, component],
                           rate_intervals[1, :, component]-rates[:, component]])
        axes[1].errorbar(q_squared*1e-12, rates[:, component], yerr=errors,
                         fmt='o', capsize=3, color=line.get_color())
    axes[1].set(xlabel=r'$q^2$ ($10^{12}$ m$^{-2}$)', ylabel=r'Rate (ps$^{-1}$)',
                title='Does one diffusion coefficient describe both components?')
    axes[1].legend(fontsize=11)
    figure.savefig(analysis.PLOT_DIR/'paper_model_diagnostics.png', dpi=180)
    plt.close(figure)
    output = analysis.PROJECT_DIR/'analysis/paper_model_comparison.json'
    output.write_text(json.dumps(results, indent=2)+'\n')
    print(json.dumps(results, indent=2))


if __name__ == '__main__':
    main()
