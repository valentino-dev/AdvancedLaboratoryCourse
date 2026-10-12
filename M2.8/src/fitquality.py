print('Fit Quality')
import matplotlib.pyplot as plt
import pandas as pd
import glob
import numpy as np
from olib import *

width = [(i+1)*0.1 for i in range(10)]
lables = [f'{width[i]:0.1f}' for i in range(10)]
width.pop(5)
width = np.array(width)
lables.pop(5)
spectra = [pd.read_csv(glob.glob(f'data/SlitWidth/{lables[i]}/*.csv')[0], delimiter=';').to_numpy().T for  i in range(9)]

p0s = list()

p0s.append([9e2, 3000, 527, 6, 0.2, 1.5])
def model_mix(x, c0, A, x0, gamma, alpha, sigma):
    return c0 + A * (alpha*gamma**2/((x-x0)**2+gamma**2)+(1-alpha)*np.exp(-(x-x0)**2/(2*sigma**2)))

p0s.append([9e2, 3000, 527, 6])
def model_lorentz(x, c0, A, x0, gamma):
    return c0 + A * (gamma**2/((x-x0)**2+gamma**2))

p0s.append([9e2, 3000, 527, 1.5])
def model_guass(x, c0, A, x0, sigma):
    return c0 + A * np.exp(-(x-x0)**2/(2*sigma**2))

models = [model_mix, model_lorentz, model_guass]




labels = ['Mixed', 'Lorentz', 'Gauss']
colors = ['red', 'green', 'blue']
red_chisqs = np.zeros((len(models), len(spectra)))
for i, spectrum in enumerate(spectra):
    plt.close('all')
    fig, ax = plt.subplots()
    for j, model in enumerate(models):
        x, y = spectrum
        mask = (x<560)*(x>500)
        x = x[mask]
        y = y[mask]
        table = Table(x, x*0, y, np.sqrt(y+1), None, 'Raman shift', 'Intensity', r'cm$^{-1}$', '(1)')
        X = np.linspace(x.min(), x.max(), 1000)
        errbar = None if j!=0 else True
        ax, (popt, pcov, red_chisq, sig) = plotData(ax, table, model=model, p0=p0s[j], show_p0=False, fit_color=colors[j], errorbar=errbar)
        red_chisqs[j, i] = red_chisq
    plt.savefig(f'images/data_point_{i}.pdf', dpi=500)


#red_chisq = red_chisqs.mean(1)
#red_chisq_err = red_chisqs.std(1)

plt.close('all')
fig, ax = plt.subplots()

for j in range(len(models)):
    x = np.arange(len(red_chisqs[j]))
    table = Table(x, x*0, red_chisqs[j], red_chisqs[j]*0, None, 'Data Point Number', r'$\chi^2/\text{ddof}$', '(1)', '(1)')
    ax = plotData(ax, table, color=colors[j], label=labels[j], errorbar=False)
plt.savefig('images/fitquality.pdf', dpi=500)


