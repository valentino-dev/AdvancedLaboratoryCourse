"""Fit PP and TG simultaneously with shared lifetimes and diffusivity.

Each channel has its own nonnegative amplitudes and unrestricted constant.
Variable projection solves these linear coefficients at each nonlinear step.
The objective gives each channel equal weight after range normalization.
Bootstrap resampling preserves correlations within a complete recorded scan.
"""

import argparse
import json

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import least_squares, nnls

import analyze as analysis


def load_channels(end_ps=None):
    channels = []
    for filename in analysis.TRACE_FILES:
        data = analysis.load_trace(filename)
        delay, scans = data[:, 0], data[:, 1:]
        _, excluded = analysis.find_anomalous_scans(delay, scans)
        keep = delay >= analysis.FIT_START_PS
        if end_ps is not None:
            keep &= delay <= end_ps
        scans = scans[keep][:, ~excluded] * 1e6
        mean = scans.mean(axis=1)
        channels.append({
            "name": filename, "time": delay[keep], "scans": scans,
            "mean": mean, "scale": np.ptp(mean),
            "q2": (analysis.nominal_period_and_wave_number(filename)[1] / 1e6)**2,
        })
    return channels


def decode(parameters, components):
    lifetimes = np.cumsum(np.exp(parameters[:components]))
    diffusion = np.exp(parameters[components])
    return lifetimes, diffusion


def project(channel, lifetimes, diffusion):
    # q2 is in 10^12 m^-2, time in ps, D in m^2/s: their numeric product is a rate.
    rates = 1 / lifetimes + diffusion * channel["q2"]
    basis = np.exp(-np.outer(channel["time"], rates))
    centered_basis = basis - basis.mean(axis=0)
    centered_signal = channel["mean"] - channel["mean"].mean()
    amplitudes, _ = nnls(centered_basis, centered_signal, maxiter=1000)
    offset = channel["mean"].mean() - basis.mean(axis=0) @ amplitudes
    prediction = offset + basis @ amplitudes
    return prediction, amplitudes, offset


def residual(parameters, channels, components, equal_channels=True):
    lifetimes, diffusion = decode(parameters, components)
    residuals = []
    for channel in channels:
        prediction, _, _ = project(channel, lifetimes, diffusion)
        weight = np.sqrt(len(prediction)) if equal_channels else 1
        residuals.append(
            (prediction - channel["mean"]) / channel["scale"] / weight
        )
    return np.concatenate(residuals)


def starts(components):
    seeds = (
        ([0.7, 8], [1.2, 20], [2, 500])
        if components == 2 else
        ([0.67, 5.2, 619], [0.7, 15, 600], [1, 4, 30], [0.5, 3, 1000])
    )
    return [
        np.log(np.r_[np.diff(np.r_[0, times]), diffusion])
        for times in seeds for diffusion in (0.001, 0.004, 0.015)
    ]


def fit_joint(channels, components=3, initial=None, equal_channels=True):
    lower = np.log([0.05] * components + [1e-8])
    upper = np.log([1e5] * components + [0.2])
    candidates = starts(components) if initial is None else initial
    best = None
    for seed in candidates:
        result = least_squares(
            residual, seed, bounds=(lower, upper),
            args=(channels, components, equal_channels),
            max_nfev=1500, ftol=1e-10, xtol=1e-10, gtol=1e-10,
        )
        if result.success and (best is None or result.cost < best.cost):
            best = result
    if best is None:
        raise RuntimeError("No joint fit converged")
    return best


def summarize(result, channels, components):
    lifetimes, diffusion = decode(result.x, components)
    traces = []
    for channel in channels:
        prediction, amplitudes, offset = project(channel, lifetimes, diffusion)
        traces.append({
            "name": channel["name"],
            "normalized_rms": float(np.sqrt(np.mean(
                ((prediction - channel["mean"]) / channel["scale"])**2))),
            "amplitudes_micro_file_units": amplitudes.tolist(),
            "offset_micro_file_units": float(offset),
        })
    return {
        "lifetimes_ps": lifetimes.tolist(), "diffusion_m2_s": float(diffusion),
        "objective": float(2 * result.cost), "channels": traces,
        "log_parameters": result.x.tolist(),
    }


def independent_reference(channels, joint):
    """Release the shared rate constraints, using the joint solution as a seed."""
    lifetimes, diffusion = decode(joint.x, 3)
    results = []
    predictions = []
    for channel in channels:
        effective_times = 1 / (1 / lifetimes + diffusion * channel["q2"])
        # q=0 here makes all three times freely fitted within this channel.
        independent = dict(channel, q2=0)
        seed = np.log(np.r_[np.diff(np.r_[0, effective_times]), 0.004])
        fitted = fit_joint([independent], initial=[seed] + starts(3))
        results.append(summarize(fitted, [independent], 3))
        predictions.append(project(independent, decode(fitted.x, 3)[0], 0)[0])
    return results, predictions


def bootstrap(channels, fitted, count, rng):
    samples = []
    failures = 0
    for iteration in range(count):
        resampled = []
        for channel in channels:
            scans = channel["scans"]
            indices = rng.integers(scans.shape[1], size=scans.shape[1])
            # Keep primary scales fixed so bootstrap measures repeat variation,
            # rather than changing the channel weighting with each resample.
            resampled.append(dict(channel, mean=scans[:, indices].mean(axis=1)))
        try:
            result = fit_joint(
                resampled, initial=[fitted.x, starts(3)[1]]
            )
            lifetimes, diffusion = decode(result.x, 3)
            samples.append(np.r_[lifetimes, diffusion])
        except (RuntimeError, ValueError):
            failures += 1
        if (iteration + 1) % 25 == 0:
            print(f"Bootstrap {iteration + 1}/{count}", flush=True)
    values = np.array(samples)
    return {
        "requested": count, "successful": len(samples), "failed": failures,
        "parameter_order": ["tau1_ps", "tau2_ps", "tau3_ps", "D_m2_s"],
        "CI95": np.percentile(values, [2.5, 97.5], axis=0).tolist(),
        "samples": values.tolist(),
    }


def plot(channels, joint_two, joint_three, independent_predictions):
    figure, axes = plt.subplots(5, 2, figsize=(11, 15), layout="constrained")
    for row, channel in enumerate(channels):
        time = channel["time"]
        mean = channel["mean"]
        scale = channel["scale"]
        origin = mean[-1]
        sem = channel["scans"].std(axis=1, ddof=1) / np.sqrt(channel["scans"].shape[1])
        predictions = [
            project(channel, *decode(joint_two.x, 2))[0],
            project(channel, *decode(joint_three.x, 3))[0],
            independent_predictions[row],
        ]
        left, right = axes[row]
        left.errorbar(time, (mean-origin)/scale, yerr=sem/scale,
                      fmt=".", color="0.35", markersize=3, elinewidth=0.6,
                      label="Mean ± repeat SE")
        for prediction, color, style, label in zip(
            predictions, ("#b06b20", "#136cb0", "#239347"), ("--", "-", ":"),
            ("Joint 2-exp", "Joint 3-exp", "Independent 3-exp"),
        ):
            left.plot(time, (prediction-origin)/scale, color=color,
                      linestyle=style, linewidth=1.4, label=label)
            right.plot(time, (mean-prediction)/scale, color=color,
                       linestyle=style, linewidth=1.1)
        right.axhline(0, color="0.3", linewidth=0.6)
        label = "PP" if row == 0 else f"TG {channel['name'].split('_')[2]} µm mask"
        left.set(title=label, ylabel="Normalized signal")
        right.set(title=f"{label}: residuals", ylabel="Residual / range")
        for axis in (left, right):
            axis.set_xlabel("Delay (ps)")
            if row == 0:
                axis.set_xscale("log")
                axis.set_xlabel("Delay (ps; logarithmic axis)")
        if row == 0:
            left.legend(fontsize=9, loc="upper right")
    figure.savefig(analysis.PLOT_DIR / "joint_pp_tg.png", dpi=180)
    figure.savefig(analysis.PLOT_DIR / "joint_pp_tg.pdf")
    plt.close(figure)

    # Report panels use larger text and separate PP from the four TG residuals.
    with plt.rc_context({
        "font.size": 14, "axes.labelsize": 14, "axes.titlesize": 15,
        "xtick.labelsize": 12, "ytick.labelsize": 12, "legend.fontsize": 11,
    }):
        figure, axes = plt.subplots(2, 1, figsize=(8, 7), layout="constrained")
        channel = channels[0]
        time, mean, scale = channel["time"], channel["mean"], channel["scale"]
        sem = channel["scans"].std(axis=1, ddof=1) / np.sqrt(channel["scans"].shape[1])
        axes[0].errorbar(time, (mean-mean[-1])/scale, yerr=sem/scale,
                        fmt=".", color="0.35", markersize=4, label="Mean ± repeat SE")
        predictions = [
            project(channel, *decode(joint_two.x, 2))[0],
            project(channel, *decode(joint_three.x, 3))[0],
            independent_predictions[0],
        ]
        for prediction, color, style, label in zip(
            predictions, ("#b06b20", "#136cb0", "#239347"), ("--", "-", ":"),
            ("Joint 2-exp", "Joint 3-exp", "Independent 3-exp"),
        ):
            axes[0].plot(time, (prediction-mean[-1])/scale, color=color,
                         linestyle=style, linewidth=1.7, label=label)
            axes[1].plot(time, (mean-prediction)/scale, color=color,
                         linestyle=style, linewidth=1.5)
        axes[0].set(ylabel="Normalized PP signal", title="PP: shared versus independent relaxation")
        axes[0].legend(ncol=2)
        axes[1].set(ylabel="Residual / range", title="PP residuals")
        axes[1].axhline(0, color="0.3", linewidth=0.7)
        for axis in axes:
            axis.set_xscale("log")
            axis.set_xlabel("Delay (ps; logarithmic axis)")
        figure.savefig(analysis.PLOT_DIR / "joint_pp_fit.pdf")
        figure.savefig(analysis.PLOT_DIR / "joint_pp_fit.png", dpi=180)
        plt.close(figure)

        figure, axes = plt.subplots(2, 2, figsize=(9, 7), layout="constrained")
        for row, (axis, channel) in enumerate(zip(axes.flat, channels[1:]), 1):
            time, mean, scale = channel["time"], channel["mean"], channel["scale"]
            for prediction, color, style, label in zip(
                [project(channel, *decode(joint_two.x, 2))[0],
                 project(channel, *decode(joint_three.x, 3))[0],
                 independent_predictions[row]],
                ("#b06b20", "#136cb0", "#239347"), ("--", "-", ":"),
                ("Joint 2-exp", "Joint 3-exp", "Independent 3-exp"),
            ):
                axis.plot(time, (mean-prediction)/scale, color=color,
                          linestyle=style, linewidth=1.3, label=label)
            axis.axhline(0, color="0.3", linewidth=0.7)
            axis.set(title=f"TG {channel['name'].split('_')[2]} µm mask",
                     xlabel="Delay (ps)", ylabel="Residual / range")
        handles, labels = axes.flat[0].get_legend_handles_labels()
        figure.legend(handles, labels, loc="outside upper center", ncol=3)
        figure.savefig(analysis.PLOT_DIR / "joint_tg_residuals.pdf")
        figure.savefig(analysis.PLOT_DIR / "joint_tg_residuals.png", dpi=180)
        plt.close(figure)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bootstrap", type=int, default=200)
    args = parser.parse_args()
    channels = load_channels()
    joint_two = fit_joint(channels, 2)
    joint_three = fit_joint(channels, 3)
    independent, independent_predictions = independent_reference(channels, joint_three)
    point_weighted = fit_joint(channels, 3, equal_channels=False)
    common_end = min(channel["time"].max() for channel in channels)
    common_channels = load_channels(common_end)
    common = fit_joint(common_channels, 3)
    results = {
        "model": "B_d + sum_j A_dj exp[-(1/tau_j + D*q_d^2)*t]",
        "weighting": "Sum of per-channel mean squared residuals divided by fixed signal range squared",
        "joint_two": summarize(joint_two, channels, 2),
        "joint_three": summarize(joint_three, channels, 3),
        "independent_three": independent,
        "point_weighted_three": summarize(point_weighted, channels, 3),
        "common_end_ps": float(common_end),
        "common_window_three": summarize(common, common_channels, 3),
        "uncertainty_scope": (
            "Whole retained-scan bootstrap; masks and primary normalization fixed. "
            "No delay/period calibration, model or selection uncertainty. "
            "RMS and objective diagnostics have no estimated uncertainties. "
            "Boundary-limited parameters are not established physical lifetimes."
        ),
    }
    print(json.dumps({key: results[key] for key in ("joint_two", "joint_three")}, indent=2), flush=True)
    plot(channels, joint_two, joint_three, independent_predictions)
    if args.bootstrap:
        rng = np.random.default_rng(20261021)
        results["bootstrap_three"] = bootstrap(channels, joint_three, args.bootstrap, rng)
        results["bootstrap_common_window"] = bootstrap(
            common_channels, common, args.bootstrap, rng
        )
    (analysis.PROJECT_DIR / "analysis/joint_pp_tg.json").write_text(
        json.dumps(results, indent=2) + "\n"
    )


if __name__ == "__main__":
    main()
