"""Temporary fixed-delay comparison using the TG filename lengths."""

import json

import matplotlib.pyplot as plt
import numpy as np

from analyze import DATA_DIR, PROJECT_DIR, PLOT_DIR, find_anomalous_scans


DELAY_PS = 5.0
LENGTHS_UM = (5.6, 6.5, 8.0, 12.0)


def main():
    """Plot actual first-repeat values and retained-repeat means at one delay."""
    rows = []
    for length in LENGTHS_UM:
        filename = f"GaAs_RT_{length:g}_X.dat"
        data = np.loadtxt(DATA_DIR / filename)
        matching_rows = np.flatnonzero(
            np.isclose(data[:, 0], DELAY_PS, rtol=0, atol=1e-9)
        )
        if len(matching_rows) != 1:
            raise ValueError(f"{filename}: expected one sample at {DELAY_PS} ps")
        row_index = matching_rows[0]
        scans = data[:, 1:]
        _, excluded = find_anomalous_scans(data[:, 0], scans)
        retained_values = scans[row_index, ~excluded]
        rows.append({
            "file": filename,
            "nominal_filename_length_um": length,
            "recorded_delay_ps": float(data[row_index, 0]),
            "first_scan_value": float(scans[row_index, 0]),
            "retained_mean": float(np.mean(retained_values)),
            "retained_mean_standard_error": float(
                np.std(retained_values, ddof=1) / np.sqrt(len(retained_values))
            ),
            "retained_scan_count": int(len(retained_values)),
        })

    # Scale only the vertical display, preserving the original file values.
    display_scale = 1e6
    first_values = [row["first_scan_value"] * display_scale for row in rows]
    means = [row["retained_mean"] * display_scale for row in rows]
    errors = [
        row["retained_mean_standard_error"] * display_scale for row in rows
    ]
    with plt.rc_context({
        "font.size": 14,
        "axes.labelsize": 15,
        "axes.titlesize": 16,
        "xtick.labelsize": 13,
        "ytick.labelsize": 13,
        "legend.fontsize": 12,
    }):
        fig, ax = plt.subplots(figsize=(9, 6), layout="constrained")
        ax.plot(
            LENGTHS_UM, first_values, "o--", color="#2864a0",
            markersize=8, linewidth=1.5,
            label="First recorded scan (line guides the eye)",
        )
        ax.errorbar(
            LENGTHS_UM, means, yerr=errors, fmt="s", color="#b84d12",
            markersize=7, capsize=5, elinewidth=1.5,
            label="Retained-scan mean ± standard error",
        )
        ax.set(
            xlabel=r"Nominal length in TG filename ($\mu$m)",
            ylabel=r"Recorded signal ($10^{-6}$ file units)",
            title="Temporary length comparison at recorded delay 5 ps",
            xticks=LENGTHS_UM,
            xlim=(5.2, 12.4),
        )
        ax.legend(loc="upper right")
        fig.savefig(PLOT_DIR / "temporary_length_snapshot.png", dpi=220)
        fig.savefig(PLOT_DIR / "temporary_length_snapshot.pdf")
        plt.close(fig)

    output = {
        "purpose": "Temporary test of a hypothesized phase proportional to filename length",
        "processing": "No interpolation, baseline subtraction, or normalization",
        "uncertainties": (
            "Mean errors are repeat-scan standard errors; individual-point, "
            "delay-calibration, and length-calibration uncertainties are unavailable."
        ),
        "rows": rows,
    }
    result_file = PROJECT_DIR / "analysis/temporary_length_snapshot.json"
    result_file.write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
