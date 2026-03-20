import numpy as np
from inspect import signature
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.differentiate import jacobian
from functools import partial
import olib


print("Random Coincidence Rates")


data = np.array(pd.read_csv("data/RandomMeasurement/Counter_Data.dat", delimiter="\t"))

t = 120.
t_err = 0.02

counts = data[:, 6:8].astype(np.float64)
counts_err = np.sqrt(counts+1)

counts = np.array([counts[0:13, :], np.flip(counts[13:26], axis=0)])
counts_err = np.array([counts_err[0:13, :], np.flip(counts_err[13:26], axis=0)])

counts = ((counts/counts_err**2).sum(axis=0)/(1/counts_err**2).sum(axis=0)).T
counts_err = (np.sqrt(1/(counts_err**-2).sum(axis=0))).T
print(counts_err)

angles = np.array(data[:, 0], dtype=np.float64)[0: 13]

count_rate = counts / t
count_rate_err = ((counts_err/t)**2 + (counts/t**2*t_err)**2)**(1/2)

def model1(x, a, b, c):
    return a+b*np.cos(x/360*2*3.141-c)

def model2(x, a, b):
    return a+b*x

def model3(x, a):
    return np.zeros(x.shape) + a

models = [model1, model1]
#p0=[[1700.,  40., 125/360*2*3.141], [0, 0]]
p0=[[1700.,  40., 125/360*2*3.141], [1700.,  40., 125/360*2*3.141]]
sigmas = list()
popts, pcovs = list(), list()

def get_sigmas(model, X, popt, pcov):
    j = np.array([jacobian(lambda p: partial(model, x)(*p), popt).df for x in X])
    sigma = np.sqrt(np.einsum('ij,jk,ik->i', j, pcov, j))
    return sigma

for i in [0, 1]:
    y = count_rate[i]
    y_err = count_rate_err[i]
    x = angles
    model = models[i]
    plt.xlabel(r"$\theta/°$")
    plt.ylabel(r"$N_"+f'{i+1}'+r"/\frac{1}{s}$")
    table = olib.Table(x, np.zeros(x.shape), y, y_err, '', r"$\theta/\degree$", r"$N_"+f'{i+1}'+r"/\frac{1}{s}$")
    table.saveAsPDF('latex/', f'rnd_coinc_det_{i}.tex')
    plt.errorbar(x, y, y_err, ls='none', marker='d', capsize=10)

    #popt, pcov = curve_fit(model, x, y, p0=p0[i], sigma=y_err, absolute_sigma=True)
    popt, pcov, sigma = olib.fit_function(model, x, y, y_err, p0[i], plot=True)
    popts.append(popt)
    pcovs.append(pcov)

    sigmas.append(sigma)
    
    
    #plt.show()
    plt.savefig(f'images/rnd_coinc_det_{i}.pdf', dpi=500, bbox_inches="tight")
    plt.clf()

X = np.linspace(min(angles), max(angles), 200)
sigmas = [get_sigmas(models[i], X, popts[i], pcovs[i]) for i in [0, 1]]
resolving_time = 0.03e-6
resolving_time_err = 0.004e-6
coincidence_rate = resolving_time * models[0](X, *popts[0]) * models[1](X, *popts[1])
coincidence_rate_err = ((coincidence_rate/resolving_time*resolving_time_err)**2+(coincidence_rate/models[0](X, *popts[0])*sigmas[0])**2+(coincidence_rate/models[1](X, *popts[1])*sigmas[1])**2)**(1/2)
Y = coincidence_rate
sigma = coincidence_rate_err
#print(coincidence_rate.std())
#print(sigma)
k=1

plt.xlabel(r"$\theta/°$")
plt.ylabel(r"$N_\text{acc}/\frac{1}{s}$")
plt.plot(X, Y)
plt.fill_between(X, Y-k*sigma, Y+k*sigma, alpha=0.3)
plt.savefig(f'images/rnd_coinc_rate.pdf', dpi=500, bbox_inches="tight")
#plt.show()
plt.clf()


#print(count_rate.shape)
#better coincidence rate
bcr = resolving_time * count_rate[0] * count_rate[1]
bcr_err = ((resolving_time_err * count_rate[0] * count_rate[1])**2+(resolving_time * count_rate_err[0] * count_rate[1])**2+(resolving_time * count_rate[0] * count_rate_err[1])**2)**(1/2)

bcr_mean = ((bcr/bcr_err**2).sum(axis=0)/(1/bcr_err**2).sum(axis=0)).T
bcr_mean_err = (np.sqrt(1/(bcr_err**-2).sum(axis=0))).T
print(bcr, bcr_err)

chi_sq = (1/bcr_err**2*(bcr-bcr_mean)**2).sum()
print("chisq:",chi_sq)


print(bcr_mean, bcr_mean_err)

direct_measurement = np.array(pd.read_csv("data/RandomMeasurement/Counter_Data.dat", delimiter="\t"))


# direct measurement count rate
dm_cr = direct_measurement[:, 8].astype(np.float64)
dm_cr_err = np.sqrt(dm_cr+1)

dm_cr = np.array([dm_cr[0:13], np.flip(dm_cr[13:26], axis=0)])
dm_cr_err = np.array([dm_cr_err[0:13], np.flip(dm_cr_err[13:26], axis=0)])

dm_cr = ((dm_cr/dm_cr_err**2).sum(axis=0)/(1/dm_cr_err**2).sum(axis=0)).T
dm_cr_err = (np.sqrt(1/(dm_cr_err**-2).sum(axis=0))).T

dm_cr_err = ((dm_cr_err/t)**2 + (dm_cr/t**2*t_err)**2)**(1/2)
dm_cr = dm_cr / t





plt.errorbar(angles, dm_cr, dm_cr_err, ls='none', capsize=10, marker='d')

x = angles
y = dm_cr
y_err = dm_cr_err
p0=[462000, 10000, 110]

model = model1



olib.fit_function(model, x, y, y_err, p0, plot=True)
plt.tight_layout()

plt.xlabel(r"$\theta/\degree$")
plt.ylabel(r"$N_\text{acc}/\frac{1}{s}$")
table = olib.Table(x, np.zeros(x.shape), y, y_err, 'Direct Measurment', r"$\theta$", r"$N_\text{acc}/\frac{1}{s}$")
table.saveAsPDF('latex/', 'rnd_dir_meas_data.tex')
plt.savefig('images/direct_measurement.pdf', dpi=500, bbox_inches="tight")


