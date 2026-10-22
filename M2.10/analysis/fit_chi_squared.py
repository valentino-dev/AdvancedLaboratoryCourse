"""Diagonal reduced-chi-squared diagnostics for the existing unweighted fits.

Mean-signal errors come from retained repeats. Rate errors are bootstrap
standard deviations, not half-widths of percentile confidence intervals.
No noise estimate is invented for the single spectrum or IAC records.
"""

import json

import numpy as np
from scipy.stats import linregress

import analyze as analysis
from paper_model_comparison import fit_three_exponentials, three_exponentials


def diagonal_chi_squared(observed, predicted, standard_errors, parameter_count):
    """Evaluate a diagnostic using diagonal errors without changing the fit."""
    observed = np.asarray(observed)
    errors = np.asarray(standard_errors)
    dof = len(observed) - parameter_count
    if dof <= 0 or np.any(~np.isfinite(errors)) or np.any(errors <= 0):
        raise ValueError("Positive errors and positive residual degrees of freedom are required")
    chi_squared = float(np.sum(((observed - predicted) / errors)**2))
    return {
        "chi_squared": chi_squared,
        "reduced_chi_squared": chi_squared / dof,
        "point_count": len(observed),
        "parameter_count": parameter_count,
        "degrees_of_freedom": dof,
    }


def fit_from(time, signal, start):
    """Use the existing fitter with a temporary start; keep scan masks fixed."""
    previous = analysis.FIT_START_PS
    try:
        analysis.FIT_START_PS = start
        return analysis.fit_trace(time, signal)
    finally:
        analysis.FIT_START_PS = previous


def sensitivity_diagnostics(channels, q2, rates, rate_errors):
    """Cover the fit-window and all-scan alternatives quoted in the report."""
    diagnostics = []
    time, scans = channels[0]
    mean = scans.mean(axis=1)
    sem = scans.std(axis=1, ddof=1)/np.sqrt(scans.shape[1])
    for start, end in ((.5, 33), (.5, 57), (.5, 100), (.5, 2000), (1, 2000), (2, 2000)):
        within = time <= end
        p, keep = fit_from(time[within], mean[within], start)
        diagnostic = diagonal_chi_squared(
            mean[within][keep], analysis.biexponential_with_offset(time[within][keep], *p),
            sem[within][keep], 5,
        )
        diagnostics.append(dict(diagnostic, label=f"PP, {start:g}--{end:g} ps"))

    generator = np.random.default_rng(20261023)
    for start in (.2, 1, 2):
        estimates = []
        for time, scans in channels[1:]:
            p, _ = fit_from(time, scans.mean(axis=1), start)
            estimates.append(1/(p[1]+p[3]))
        samples = []
        for _ in range(analysis.BOOTSTRAP_REPLICATES):
            iteration = []
            for time, scans in channels[1:]:
                columns = generator.integers(0, scans.shape[1], scans.shape[1])
                p, _ = fit_from(time, scans[:, columns].mean(axis=1), start)
                iteration.append(1/(p[1]+p[3]))
            samples.append(iteration)
        errors = np.std(samples, axis=0, ddof=1)
        regression = linregress(q2[1:], estimates)
        diagnostic = diagonal_chi_squared(
            estimates, regression.intercept+regression.slope*q2[1:], errors, 2,
        )
        diagnostics.append(dict(diagnostic, label=f"TG slow-rate line, start {start:g} ps"))
        print(f"Completed chi-squared window check: TG start {start:g} ps", flush=True)

    # The other scan-selection alternative quoted in the report retains both
    # previously excluded mask-8 scans; other datasets remain unchanged.
    data8 = analysis.load_trace("GaAs_RT_8_X.dat")
    time8, scans8 = data8[:, 0], data8[:, 1:]
    p, _ = analysis.fit_trace(time8, scans8.mean(axis=1))
    all_rates = rates[:, 1].copy()
    all_rates[3] = 1/(p[1]+p[3])
    all_errors = rate_errors[:, 1].copy()
    samples = []
    for _ in range(analysis.BOOTSTRAP_REPLICATES):
        columns = generator.integers(0, scans8.shape[1], scans8.shape[1])
        p, _ = analysis.fit_trace(time8, scans8[:, columns].mean(axis=1))
        samples.append(1/(p[1]+p[3]))
    all_errors[3] = np.std(samples, ddof=1)
    regression = linregress(q2[1:], all_rates[1:])
    diagnostic = diagonal_chi_squared(
        all_rates[1:], regression.intercept+regression.slope*q2[1:], all_errors[1:], 2,
    )
    diagnostics.append(dict(diagnostic, label="TG slow-rate line, all scans retained"))
    return diagnostics


def calculate():
    channels, fits, transients = [], [], {}
    for filename in analysis.TRACE_FILES:
        data = analysis.load_trace(filename)
        time, scans = data[:, 0], data[:, 1:]
        _, excluded = analysis.find_anomalous_scans(time, scans)
        scans = scans[:, ~excluded]
        mean = scans.mean(axis=1)
        sem = scans.std(axis=1, ddof=1) / np.sqrt(scans.shape[1])
        p, keep = analysis.fit_trace(time, mean)
        prediction = analysis.biexponential_with_offset(time[keep], *p)
        diagnostic = diagonal_chi_squared(mean[keep], prediction, sem[keep], 5)
        diagnostic["error_source"] = "Pointwise standard error of retained-scan mean"
        transients[filename] = diagnostic
        fits.append(p)
        channels.append((time, scans))

    pp_time, pp_scans = channels[0]
    p3, _, _, keep = fit_three_exponentials(pp_time, pp_scans.mean(axis=1))
    pp_sem = pp_scans.std(axis=1, ddof=1) / np.sqrt(pp_scans.shape[1])
    pp_three = diagonal_chi_squared(
        pp_scans.mean(axis=1)[keep],
        three_exponentials(pp_time[keep], *p3)*1e-6,
        pp_sem[keep], 7,
    )

    data8 = analysis.load_trace("GaAs_RT_8_X.dat")
    retained8 = np.ones(data8.shape[1]-1, dtype=bool)
    retained8[9] = False
    sensitivity_scans = data8[:, 1:][:, retained8]
    p8, _ = analysis.fit_trace(data8[:, 0], sensitivity_scans.mean(axis=1))

    # Repeat the primary bootstrap sampling sequence. A separate generator for
    # the scan-4 sensitivity prevents that extra calculation changing it.
    generator = np.random.default_rng(20261019)
    sensitivity_generator = np.random.default_rng(20261022)
    bootstrap_rates, sensitivity_samples = [], []
    for _ in range(analysis.BOOTSTRAP_REPLICATES):
        iteration = []
        for time, scans in channels:
            columns = generator.integers(0, scans.shape[1], scans.shape[1])
            p, _ = analysis.fit_trace(time, scans[:, columns].mean(axis=1))
            iteration.append([1/p[1], 1/(p[1]+p[3])])
        bootstrap_rates.append(iteration)
        columns = sensitivity_generator.integers(
            0, sensitivity_scans.shape[1], sensitivity_scans.shape[1]
        )
        p, _ = analysis.fit_trace(data8[:, 0], sensitivity_scans[:, columns].mean(axis=1))
        sensitivity_samples.append(1/(p[1]+p[3]))

    rate_errors = np.std(bootstrap_rates, axis=0, ddof=1)
    rates = np.array([[1/p[1], 1/(p[1]+p[3])] for p in fits])
    q2 = np.array([analysis.nominal_period_and_wave_number(name)[1]**2
                   for name in analysis.TRACE_FILES])
    sensitivity_rates = rates[:, 1].copy()
    sensitivity_rates[3] = 1/(p8[1]+p8[3])
    sensitivity_errors = rate_errors[:, 1].copy()
    sensitivity_errors[3] = np.std(sensitivity_samples, ddof=1)
    specifications = [
        ("slow_tg", q2[1:], rates[1:, 1], rate_errors[1:, 1]),
        ("slow_with_pp", q2, rates[:, 1], rate_errors[:, 1]),
        ("slow_keep_scan4", q2[1:], sensitivity_rates[1:], sensitivity_errors[1:]),
        ("fast_tg", q2[1:], rates[1:, 0], rate_errors[1:, 0]),
    ]
    regressions = {}
    for name, x, y, errors in specifications:
        fit = linregress(x, y)
        regressions[name] = diagonal_chi_squared(y, fit.intercept+fit.slope*x, errors, 2)
        regressions[name]["standard_errors_per_ps"] = errors.tolist()
        regressions[name]["error_source"] = "Whole-scan bootstrap standard deviation of fitted rates"

    result = {
        "transients": transients, "pp_three": pp_three, "regressions": regressions,
        "sensitivity_fits": sensitivity_diagnostics(channels, q2, rates, rate_errors),
        "unavailable": {
            "spectrum": "No pointwise measurement uncertainties or independent noise model",
            "iac": "No pointwise measurement uncertainties or independent noise model",
            "fringe_derivative": "Not independently fitted; derived from the IAC fit",
        },
        "interpretation": (
            "Descriptive diagonal diagnostics for existing unweighted fits. "
            "Correlations between delay points are omitted and error scales are estimated. "
            "N-p is a nominal denominator, not a calibrated reference distribution; "
            "no formal chi-squared goodness-of-fit probability is inferred. "
            "Diagnostic uncertainties are not estimated. Delays and spacings are exact inputs."
        ),
    }
    output = analysis.PROJECT_DIR / "analysis/fit_chi_squared.json"
    output.write_text(json.dumps(result, indent=2)+"\n")
    return result


if __name__ == "__main__":
    print(json.dumps(calculate(), indent=2))
