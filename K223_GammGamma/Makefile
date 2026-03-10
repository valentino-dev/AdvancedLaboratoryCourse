# Python interpreter
PYTHON := python3

# Main script
SCRIPT := src/prompt_curve.py

.PHONY: prompt_curve spectrum threshold_scan angular correlation

prompt_curve:
	$(PYTHON) $(SCRIPT)
spectrum:
	$(PYTHON) src/pseudo/spectrum.py
threshold_scan:
	$(PYTHON) src/pseudo/threshold_scan.py
angular_correlation:
	$(PYTHON) src/pseudo/angular_correlation.py
