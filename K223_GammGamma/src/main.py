import numpy as np
import olib
from inspect import signature
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.differentiate import jacobian
from functools import partial
from scipy.special import legendre

print("Angular Correlation Function")

#data = np.array(pd.read_csv("data/MainMeasurement/Config.ini", delimiter="\t"))
data = np.array(pd.read_csv("data/MainMeasurement/Counter_Data.dat", delimiter="\t"))

t = 12*60.
t_err = 0.02

angles = np.array(data[:, 0], dtype=np.float64)[0: 13]

#data = data[:, 11:13].astype(np.float64).T
data = data[:, 7:9].astype(np.float64).T

def trans(counts):
    return np.array([counts[:, i*len(angles):(i+1)*len(angles)] if i%2==0 else np.flip(counts[:, i*len(angles):(i+1)*len(angles)], axis=1) for i in range(counts.shape[1]//len(angles))]).transpose((1, 0, 2))

counts = trans(data)
counts_err = np.sqrt(counts+1)

count_rate = counts/t
count_rate_err = ((counts_err/t)**2+(counts/t**2*t_err)**2)**(1/2)

count_rate = np.average(count_rate, axis=1, weights=count_rate_err**-2)
count_rate_err = np.sqrt(1/np.sum(count_rate_err**-2, axis=1))


def model(x, N, A22, A44):
    return N*(1+A22 * legendre(2)(np.cos(x/360*2*3.141))+ A44 *legendre(4)(np.cos(x/360*2*3.141)))

def distance_sq(r_d, r_s, theta_d, theta_s):
    return r_d**2-2*r_s*r_d*np.cos(theta_d - theta_s)+r_s**2

def corr_model(theta_d, r_s, theta_s):
    r_d = 10e-3/2
    theta_d = theta_d/360*2*3.141
    theta_s = theta_s/360*2*3.141
    return distance_sq(r_d, r_s, theta_d, theta_s)/distance_sq(r_d, r_s, 180/360*2*3.141, theta_s)

def model1(x, a, b, c):
    return a+b*np.cos(x/360*2*3.141-c)

x = angles
y = count_rate[1]*corr_model(angles, 64.9e-6, 116.5) - model1(angles, 0.153, -0.025, 110.87)
y_err = count_rate_err[1]
plt.errorbar(x, y, y_err, ls='none', marker='d', capsize=10)

p0=[4, 0.1020, 0.0091]
plt.xlabel(r"$\theta/°$")
plt.ylabel(r"Corrected Coincidence Rate /$\frac{1}{s}$")

popt, pcov, _ = olib.fit_function(model, x, y, y_err, p0, plot=True, show_p0=False)
plt.savefig(f'images/main.pdf', dpi=500, bbox_inches="tight")

print(*olib.roundToError(*(popt[1:]/np.array([0.9057**2, 0.7113**2]))))
print(*olib.roundToError(*(np.sqrt(np.diag(pcov))[1:]/np.array([0.9057**2, 0.7113**2]))))

plt.clf()
plt.xlabel(r"$\theta/°$")
plt.ylabel(r"Corrected Coincidence Rate /$\frac{1}{s}$")
plt.errorbar(x, y, y_err, ls='none', marker='d', capsize=10)
F_22_1 = np.array([0.707, -0.598, 0.071, -0.177, -0.171])
F_44_1 = np.array([0, -1.069, 0, 0, -0.008])
F_22_2 = np.array([0.707, -0.598, 0.707, 0.707, -0.598])
F_44_2 = np.array([0, -1.069, 0, 0, -1.069])
A_22 = F_22_1*F_22_2*0.9057**2
A_44 = F_44_1*F_44_2*0.7113**2
lables = ['0(1)1(1)0', '0(2)2(2)0', '2(1)1(1)0', '4(3)1(1)0', '4(2)2(2)0']

for i in range(A_22.shape[0]):
    X = np.linspace(min(angles), max(angles), 200)
    Y = model(X, popt[0], A_22[i], A_44[i])
    olib.fit_function(lambda x, N: model(x, N, A_22[i], A_44[i]), angles, y, y_err, [popt[0]], plot=True, label=lables[i])

plt.legend()
plt.savefig(f'images/main_other.pdf', dpi=500, bbox_inches="tight")
