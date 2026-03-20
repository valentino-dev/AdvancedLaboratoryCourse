# Python interpreter
PYTHON := python3

# Main script
SCRIPT := src/prompt_curve.py

.PHONY: prompt_curve spectrum threshold_scan angular correlation sca_calib rnd_coinc_rate stability misal main

prompt_curve:
	$(PYTHON) $(SCRIPT)
spectrum:
	$(PYTHON) src/pseudo/spectrum.py
threshold_scan:
	$(PYTHON) src/pseudo/threshold_scan.py
angular_correlation:
	$(PYTHON) src/pseudo/angular_correlation.py

sca_calib:
	python3 src/sca_calib.py

rnd_coinc_rate:
	python3 src/rnd_coinc_rate.py

stab:
	python3 src/stability.py

misal:
	python3 src/misal.py

main:
	python3 src/main.py
