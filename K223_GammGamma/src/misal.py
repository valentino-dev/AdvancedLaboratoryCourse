import numpy as np
import olib
from inspect import signature
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.differentiate import jacobian
from functools import partial

print("Missalignment check")

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

count_rate = np.average(count_rate, axis=1, weights=count_rate_err**-2)
count_rate_err = np.sqrt(1/np.sum(count_rate_err**-2, axis=1))

count_rate = count_rate[0]
count_rate_err = count_rate_err[0]

ratio = count_rate[6]/count_rate
ratio_err = ((count_rate[6]/count_rate**2 * count_rate_err)**2+(count_rate_err[6]/count_rate)**2)**(1/2)
x = angles
y = ratio
y_err = ratio_err

def distance_sq(r_d, r_s, theta_d, theta_s):
    return r_d**2-2*r_s*r_d*np.cos(theta_d - theta_s)+r_s**2

def model(theta_d, r_s, theta_s):
    r_d = 10e-3/2
    theta_d = theta_d/360*2*3.141
    theta_s = theta_s/360*2*3.141
    return distance_sq(r_d, r_s, theta_d, theta_s)/distance_sq(r_d, r_s, 180/360*2*3.141, theta_s)

plt.errorbar(x, y, y_err, ls='none', marker='d', capsize=10)
p0 = [1.4e-4, 120]
_ = olib.fit_function(model, x, y, y_err, p0, plot=True)
plt.xlabel(r"$\theta/°$")
plt.ylabel(r"Ratio $\frac{N(180°)}{N(\theta)}$")
plt.savefig(f'images/misal.pdf', dpi=500, bbox_inches="tight")
