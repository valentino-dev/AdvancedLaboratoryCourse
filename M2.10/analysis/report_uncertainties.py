"""Propagate conditional uncertainties for quantities quoted in the report.

Run after analyze.py and review_checks.py. Calibration and model uncertainties
are unavailable; they must not be represented by these statistical errors.
"""
import json
from pathlib import Path

import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import linregress

import analyze as analysis


def standard_error(gradient, covariance):
    gradient = np.asarray(gradient)
    return float(np.sqrt(gradient @ covariance @ gradient))


def main():
    primary = json.loads((analysis.PROJECT_DIR / 'analysis/results.json').read_text())
    checks = json.loads((analysis.PROJECT_DIR / 'analysis/review_checks.json').read_text())
    x, y = np.loadtxt(analysis.DATA_DIR / 'spektrum.csv', delimiter=',').T
    p, covariance = curve_fit(analysis.gaussian_with_offset, x, y,
                              p0=(y.max()-np.median(y), x[np.argmax(y)], 5, np.median(y)),
                              maxfev=20000)
    amplitude, wavelength, sigma, background = p
    width = analysis.GAUSSIAN_FWHM * abs(sigma)
    frequency_thz = analysis.SPEED_OF_LIGHT * width / wavelength**2 * 1e-3
    frequency_gradient = [0, -2*frequency_thz/wavelength, frequency_thz/sigma, 0]
    frequency_error = standard_error(frequency_gradient, covariance)
    limit_fs = analysis.GAUSSIAN_TIME_BANDWIDTH / frequency_thz * 1e3
    results = {
        'wavelength_nm': [wavelength, np.sqrt(covariance[1, 1])],
        'peak_level': [background+amplitude, standard_error([1,0,0,1],covariance)],
        'half_maximum_level': [background+amplitude/2, standard_error([.5,0,0,1],covariance)],
        'bandwidth_nm': [width, analysis.GAUSSIAN_FWHM*np.sqrt(covariance[2,2])],
        'frequency_thz': [frequency_thz, frequency_error],
        'transform_limit_fs': [limit_fs, limit_fs*frequency_error/frequency_thz],
    }
    selected = {}
    for group, name in primary['selected_autocorrelation'].items():
        fit = next(item for item in primary['autocorrelation'][group] if item['filename']==name)
        selected[group] = fit
        pulse = fit['pulse_width_fs']
        conditional_error = fit['pulse_stderr_fs']
        # The pulse width is proportional to lambda; spectrum and IAC records
        # are treated as independent, but derived spectral quantities are not.
        pulse_error = np.hypot(conditional_error, pulse*np.sqrt(covariance[1,1])/wavelength)
        product = pulse*frequency_thz*1e-3
        # TBP ∝ fitted IAC width in cycles * spectral sigma / lambda.
        # This accounts for the shared wavelength in pulse width and bandwidth.
        spectral_gradient = [0, -product/wavelength, product/sigma, 0]
        product_error = np.hypot(product*conditional_error/pulse,
                                 standard_error(spectral_gradient,covariance))
        results[group+'_pulse_fs'] = [pulse, pulse_error]
        results[group+'_time_bandwidth'] = [product, product_error]
    before, after = selected['vor'], selected['nach']
    ratio = after['pulse_width_fs']/before['pulse_width_fs']
    ratio_error = ratio*np.hypot(after['pulse_stderr_fs']/after['pulse_width_fs'],
                                 before['pulse_stderr_fs']/before['pulse_width_fs'])
    results['compression_percent'] = [100*(1-ratio), 100*ratio_error]

    q_squared = np.asarray(primary['delay']['q_squared'])
    rates = np.asarray(primary['delay']['rates'])
    regression = linregress(q_squared[1:], rates[1:])
    results['tg_intercept_per_ps'] = [regression.intercept, regression.intercept_stderr]
    results['nominal_recombination_ps'] = [1/regression.intercept,
                                          regression.intercept_stderr/regression.intercept**2]
    # T is an assumed fixed value, with no measured uncertainty available.
    mobility_factor = 1.602176634e-19/(1.380649e-23*300)*1e4
    results['mobility_cm2_Vs'] = [regression.slope*1e12*mobility_factor,
                                 regression.stderr*1e12*mobility_factor]
    results['mobility_bootstrap_ci'] = (np.array(primary['delay']['tg_bootstrap_ci'])*mobility_factor).tolist()
    ci = np.array(primary['delay']['lifetime_bootstrap_ci'])
    results['pp_slow_rate_ci'] = (1/ci[::-1,0,1]).tolist()

    results['spectrum_windows'] = []
    for half_window in (50,30,20):
        keep = abs(x-789)<half_window
        params, cov = curve_fit(analysis.gaussian_with_offset,x[keep],y[keep],p0=(11000,789,5,2500))
        results['spectrum_windows'].append([half_window,analysis.GAUSSIAN_FWHM*abs(params[2]),
                                            analysis.GAUSSIAN_FWHM*np.sqrt(cov[2,2])])
    results['diffusion_start_sensitivity'] = []
    for check in checks['delay_start']:
        lifetimes = np.array(check['lifetimes_ps'])
        fit = linregress(q_squared[1:],1/lifetimes[1:,1])
        results['diffusion_start_sensitivity'].append([check['start_ps'],fit.slope*1e12,fit.stderr*1e12])
    for label, retain_scan4, retain_scan10 in [('retain_scan4',True,False),('all_scans',True,True)]:
        data=analysis.load_trace('GaAs_RT_8_X.dat')
        keep=np.ones(data.shape[1]-1,dtype=bool)
        keep[3],keep[9]=retain_scan4,retain_scan10
        parameters,_=analysis.fit_trace(data[:,0],data[:,1:][:,keep].mean(axis=1))
        changed_rates=rates.copy()
        changed_rates[3]=1/(parameters[1]+parameters[3])
        fit=linregress(q_squared[1:],changed_rates[1:])
        results[label+'_diffusion']=[fit.slope*1e12,fit.stderr*1e12]
    data=analysis.load_trace('GaAs_RT_PP_X.dat')
    time, mean=data[:,0],data[:,1:].mean(axis=1)
    results['pp_window_covariance_errors']=[]
    original_start=analysis.FIT_START_PS
    try:
        for start,end in [(0.5,33),(0.5,57),(0.5,100),(0.5,2000),(1,2000),(2,2000)]:
            analysis.FIT_START_PS=start
            keep=time<=end
            params,_,cov=analysis.fit_trace(time[keep],mean[keep],return_covariance=True)
            results['pp_window_covariance_errors'].append({
                'start_ps':start,'end_ps':end,
                'fast_ps':[params[1],np.sqrt(cov[1,1])],
                'slow_ps':[params[1]+params[3],standard_error([0,1,0,1,0],cov)]})
    finally:
        analysis.FIT_START_PS=original_start
    results['iac_window_endpoints']=[]
    for group,name in primary['selected_autocorrelation'].items():
        for window in (2.0,3.2):
            fit=analysis.fit_autocorrelation(name,wavelength,window_ms=window)
            results['iac_window_endpoints'].append({
                'name':name,'half_window_ms':window,
                'pulse_fs':[fit['pulse_width_fs'],fit['pulse_stderr_fs']]})
    results['constant_scan_speed']=[]
    for group,name in primary['selected_autocorrelation'].items():
        fit=analysis.fit_autocorrelation(name,wavelength,nonlinear_scan=False)
        results['constant_scan_speed'].append({
            'name':name,'pulse_fs':[fit['pulse_width_fs'],fit['pulse_stderr_fs']]})
    results['conventions'] = {
        'pairs': '[estimate, conditional one-standard-error uncertainty]',
        'ci': '95% retained-whole-scan bootstrap percentile intervals',
        'exclusions': 'Supplied PP/TG delays and grating spacings are exact inputs. Wavelength calibration, model, selection and assumed-temperature uncertainty are not included.',
        'dependence': 'Full spectral covariance retained; common wavelength cancels in compression ratio.'
    }
    (analysis.PROJECT_DIR/'analysis/uncertainties.json').write_text(json.dumps(results,indent=2)+'\n')
    print(json.dumps(results,indent=2))


if __name__=='__main__':
    main()
