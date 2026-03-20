import numpy as np
from inspect import signature
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.differentiate import jacobian
from functools import partial
import olib

print("Stability check")

data = np.array(pd.read_csv("data/MainMeasurement/Counter_Data.dat", delimiter="\t"))

t = 12*60.
t_err = 0.02

angles = np.array(data[:, 0], dtype=np.float64)[0: 13]

data = data[:, 6:8].astype(np.float64).T

def trans(counts):
    return np.array([counts[:, i*len(angles):(i+1)*len(angles)] if i%2==0 else np.flip(counts[:, i*len(angles):(i+1)*len(angles)], axis=1) for i in range(counts.shape[1]//len(angles))]).transpose((1, 0, 2))

counts = trans(data)
counts_err = np.sqrt(counts+1)

count_rate = counts/t
count_rate_err = ((counts_err/t)**2+(counts/t**2*t_err)**2)**(1/2)

mean = np.average(count_rate, axis=1, weights=count_rate_err**-2)
mean_err = np.sqrt(1/np.sum(count_rate_err**-2, axis=1))
std = count_rate.std(axis=1)


for i in [0, 1]:
    plt.clf()

    relative_std = std[i]/mean[i]
    relative_std_err = std[i]/mean[i]**2*mean_err[i]

    relative_uncertainty = mean_err[i]/mean[i]
    plt.scatter(angles, relative_std, marker='d', label='relative std')
    plt.scatter(angles, relative_uncertainty, marker='d', label='relative uncertainty')
    x = angles
    y = relative_std
    y_err = np.zeros(y.shape)
    table = olib.Table(x, np.zeros(x.shape), y, y_err, 'Direct Measurment', r"$\theta$", r"$N_\text{acc}/\frac{1}{s}$")
    table.saveAsPDF('latex/', f"stab_count_det_err{i+1}.tex")
    y = relative_uncertainty
    y_err = np.zeros(y.shape)
    table = olib.Table(x, np.zeros(x.shape), y, y_err, 'Direct Measurment', r"$\theta$", r"$N_\text{acc}/\frac{1}{s}$")
    table.saveAsPDF('latex/', f"stab_count_det_unc{i+1}.tex")
    plt.xlabel(r"$\theta/°$")
    plt.ylabel(r"relative standard deviation or relative unceratinty")
    plt.legend()
    plt.yscale('log')
    plt.savefig(f"images/stab_count_det{i+1}.pdf", dpi=500, bbox_inches="tight")

    plt.clf()
    plt.xlabel(r"time $t/\text{min.}$")
    plt.ylabel(r'count rate $N/\frac{1}{\text{s}}$')
    lss=['-']*13
    lss[10]='--'
    lss[11]='--'
    lss[12]='--'
    for k in range(13):
        times = [j*12*13+k*12 if j%2==0 else (j+1)*12*13-(k+1)*12 for j in range(6)]
        plt.plot(times, count_rate[i, :, k]/count_rate_err[i, 0, k], label=f'$\\theta={k*15+90}°$', ls=lss[k])

    plt.legend()
    plt.savefig(f"images/stab_count_det_behaviour_{i+1}.pdf", dpi=500, bbox_inches="tight")

data = np.array(pd.read_csv("data/MainMeasurement/Counter_Data.dat", delimiter="\t"))

data = data[:, 8:10].astype(np.float64).T
def trans(counts):
    return np.array([counts[:, i*len(angles):(i+1)*len(angles)] if i%2==0 else np.flip(counts[:, i*len(angles):(i+1)*len(angles)], axis=1) for i in range(counts.shape[1]//len(angles))]).transpose((1, 0, 2))

counts = trans(data)
counts_err = np.sqrt(counts+1)
#print(counts_err/counts)

count_rate = counts/t
count_rate_err = ((counts_err/t)**2+(counts/t**2*t_err)**2)**(1/2)

mean = np.average(count_rate, axis=1, weights=count_rate_err**-2)
mean_err = np.sqrt(1/np.sum(count_rate_err**-2, axis=1))
std = count_rate.std(axis=1)

for i in [0, 1]:
    plt.clf()

    relative_std = std[i]/mean[i]
    relative_std_err = std[i]/mean[i]**2*mean_err[i]

    relative_uncertainty = mean_err[i]/mean[i]
    plt.scatter(angles, relative_std, marker='d', label='relative std')
    #plt.errorbar(angles, relative_std, relative_std_err, marker='d', ls='none', capsize=10)
    plt.scatter(angles, relative_uncertainty, marker='d', label='relative uncertainty')
    x = angles
    y = relative_std
    y_err = np.zeros(y.shape)
    table = olib.Table(x, np.zeros(x.shape), y, y_err, 'Direct Measurment', r"$\theta$", r"$N_\text{acc}/\frac{1}{s}$")
    table.saveAsPDF('latex/', f"stab_count_circuit_err{i+1}.tex")
    y = relative_uncertainty
    y_err = np.zeros(y.shape)
    table = olib.Table(x, np.zeros(x.shape), y, y_err, 'Direct Measurment', r"$\theta$", r"$N_\text{acc}/\frac{1}{s}$")
    table.saveAsPDF('latex/', f"stab_count_circuit_unc{i+1}.tex")
    plt.xlabel(r"$\theta/°$")
    plt.ylabel(r"relative standard deviation or relative unceratinty")
    plt.legend()
    plt.yscale('log')
    plt.savefig(f"images/stab_count_circuit{i+1}.pdf", dpi=500, bbox_inches="tight")

    plt.clf()
    plt.xlabel(r"time $t/\text{min.}$")
    plt.ylabel(r'count rate $N/\frac{1}{\text{s}}$')
    lss=['-']*13
    lss[10]='--'
    lss[11]='--'
    lss[12]='--'
    for k in range(13):
        times = [j*12*13+k*12 if j%2==0 else (j+1)*12*13-(k+1)*12 for j in range(6)]
        plt.plot(times, count_rate[i, :, k]/count_rate_err[i, 0, k], label=f'$\\theta={k*15+90}°$', ls=lss[k])

    plt.legend()
    plt.savefig(f"images/stab_count_circuit_behaviour_{i+1}.pdf", dpi=500, bbox_inches="tight")
