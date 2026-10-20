# Analysis review — 19 October 2026

## Corrected issues

1. **Autocorrelation width conversion.** The earlier fit assigned the same Gaussian width to the DC, fundamental and second-harmonic terms. Direct evaluation of the tutorial integrals gives different widths. The fit now preserves those relationships and fits the pulse-intensity standard deviation in optical cycles directly.
2. **Nonuniform delay scan.** Observed fringe spacings change across the pulse overlap. A cubic mapping from recorded time to optical cycles replaces the constant-speed approximation. Sine/cosine coefficients avoid artificial phase-boundary constraints. Every raw sample in the selected window enters the fit; no Hilbert transform is used.
3. **Tutorial expansion.** Printed Eq. (4.2) omits relative numerical factors present in the expansion of Eq. (4.1). The report retains the printed equation, derives the corrected coefficients, and checks the ideal 8:1 peak/background ratio.
4. **Envelope interpretation.** The requested square-root overlay is explicitly identified as the square root of an autocorrelation modulation bound. It is neither the raw-signal upper envelope nor the laser electric-field envelope. Width extraction uses the raw fit, not this display curve.
5. **Exclusion reasoning.** `vor2` remains excluded because substantial oscillations persist at its record boundary. Its different sampling interval alone is not grounds for exclusion. TG exclusions are statistical choices, not demonstrated instrument failures; both exclusion sensitivities are reported.
6. **Physical interpretation.** Corrected the reversed rate/lifetime comparison, distinguished thermalization from cooling, and derived the conditional relation between cooling time and electron–phonon coupling. No unsupported numerical volumetric coupling coefficient is claimed.
7. **Uncertainty and plotting.** Added fit-window and spectral-window sensitivity calculations; stated that transient fits are unweighted. Fit lines now match data colors. TG residuals no longer share the 2000 ps PP axis. The spectrum has a peak detail panel; rate plots show repeat-scan intervals and the actual borderline-scan sensitivity line. Fixed zero-MAD handling and the excluded-scan plot's single-scan-axis edge case.

## Numerical consequences

- Selected IAC traces: `vor3`, **112.894 fs**, and `nach2`, **93.902 fs**. Normalized RMS: 0.02227 and 0.01342. These are effective Gaussian-model estimates, not independently calibrated absolute pulse measurements.
- Primary spectrum: 788.656 nm center, 10.916 nm FWHM. Gaussian transform limit: approximately 83.9 fs. Time–bandwidth products: approximately 0.594 and 0.494.
- IAC half-windows of 2.0–3.2 ms give 111.7–114.1 fs before and 92.7–94.3 fs after compression. These variations exceed formal covariance errors.
- The TG lifetimes, quarter-period conversion and diffusion conversion were confirmed. Primary diffusivity remains **0.0037134 m²/s** (37.13 cm²/s).
- Retaining all TG scans gives 0.0037778 m²/s; the main transport conclusion does not depend strongly on the exclusions.
- PP slow lifetime changes from 2.93 ps to 7.88 ps as the fit end changes from 33 ps to 2000 ps. A common microscopic assignment is not supported by the small bootstrap intervals alone.

## Checks performed

- Recomputed the full analysis, including 300 whole-scan bootstrap replicates.
- Numerically integrated the defining fourth-power field expression and confirmed the derived IAC formula at multiple delays.
- Tried five starting lifetime pairs for each primary transient fit; they converged to the same reported solution within numerical precision.
- Checked identical and zero-valued repeat arrays: finite screening scores and no automatic exclusions.
- Evaluated spectral fitting windows, IAC fitting windows, constant versus nonlinear scan calibration, PP/TG fitting windows, and inclusion of flagged TG scans.
- Compiled the 23-page report with resolved references and no overfull boxes; inspected rendered equations, tables and figures. Both published workspace PDF copies have identical hashes.
- Verified the cited numerical literature comparisons against the primary papers' abstracts; linked sources remain in the report bibliography.

## Reproduce

From the project root:

```bash
OPENBLAS_NUM_THREADS=1 MPLBACKEND=Agg python analysis/analyze.py > analysis/results.txt
OPENBLAS_NUM_THREADS=1 MPLBACKEND=Agg python analysis/review_checks.py
OPENBLAS_NUM_THREADS=1 MPLBACKEND=Agg python analysis/report_uncertainties.py
MPLBACKEND=Agg python analysis/plot_vor_nach.py
```

`analyze.py` writes plots, `analysis/results.json`, and the generated LaTeX candidate table `latex/autocorrelation_results.tex`. `review_checks.py` writes `analysis/review_checks.json` and the fringe-calibration diagnostic. The unused `olib.py` is not part of the active calculation path.

Compile from `latex/` with `pdflatex -interaction=nonstopmode -halt-on-error main.tex` three times to resolve references and the table of contents. The reviewed build was made in `latex/build-review/`; its PDF is copied to `latex/main.pdf` and the project-root `main.pdf`.

## Remaining limitations

Unknown detector filtering, pulse chirp, excitation density, actual sample temperature and independent geometrical/time calibration remain experimental limitations. The cubic scan mapping and Gaussian overlap model improve agreement with the data but do not uniquely reconstruct the pulse's spectral phase. Bootstrap and covariance intervals remain conditional on their models and acquisition assumptions.

## Uncertainty reporting update

The report now pairs physical estimates with conditional standard errors or explicitly identified confidence intervals, including repeated values in the abstract and conclusion. Propagated errors for spectral quantities, time–bandwidth products, compression, the TG intercept and mobility are recorded in `uncertainties.json`. The shared wavelength dependence is included. Diagnostic statistics, nominal calibrations and literature values without a supplied uncertainty are explicitly distinguished. Alternative PP-window errors are conditional fit covariance errors, not the repeat-scan bootstrap intervals used for the primary transient results.

## Additional paper: Chen et al. (arXiv:1805.01561v1)

The newly supplied paper motivates channel-specific signal models and a paired-bootstrap comparison of fast/slow diffusion slopes. `paper_model_comparison.py` writes `paper_model_comparison.json` and `plots/paper_model_diagnostics.png`. The report adds the phase-separation derivation, the shared-slope check, and a three-exponential PP sensitivity fit. The alternative PP lifetimes have conditional covariance errors; they do not identify the slow contribution as recombination. The primary fit remains explicitly a baseline model.
