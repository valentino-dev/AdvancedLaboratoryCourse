"""Plot the six auxiliary before/after traces in Data TG 4_bis."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


PROJECT_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = PROJECT_DIR / "Data TG 4_bis"
OUTPUT_FILE = PROJECT_DIR / "plots" / "vor_nach.png"


def plot_group(axis, names, title):
    """Plot each two-column time/signal file in one panel."""
    for name in names:
        time_seconds, signal = np.loadtxt(DATA_DIR / name, delimiter=",").T
        axis.plot(time_seconds * 1e3, signal, linewidth=0.7, alpha=0.8, label=name)

    axis.set_title(title)
    axis.set_xlabel("Recorded time (ms)")
    axis.set_ylabel("Signal (recorded units)")
    axis.grid(True, alpha=0.3)
    axis.legend()


def main():
    figure, axes = plt.subplots(2, 1, figsize=(10, 7), constrained_layout=True)
    plot_group(axes[0], ("vor1", "vor2", "vor3"), "Before (vor)")
    plot_group(axes[1], ("nach1", "nach2", "nach3"), "After (nach)")
    figure.suptitle("Recorded vor and nach traces")
    OUTPUT_FILE.parent.mkdir(exist_ok=True)
    figure.savefig(OUTPUT_FILE, dpi=180)
    print(f"Saved plot to {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
