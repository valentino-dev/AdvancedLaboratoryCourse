import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import least_squares

print("Spectrum")


data = np.loadtxt('data/pseudo/spectrum.dat', delimiter='\t').T

def gaussian_fct(sig, mu, x):
    return 1/(sig*(2*3.141)**(1/2))*np.exp(-(x-mu)**2/(2*sig**2))
def fct(A, sig, mu, m, b, x):
    return A * (sig*(2*3.141)**(1/2))* gaussian_fct(sig, mu, x) + m*x +b

p0_1 = [6000, 10, 127, 0, 0]
p0_2 = [5000, 7, 146, 0, 0]

def residuals(params, x, y):
    return fct(*params, *[x]) - y

l_lim_1 = data[0]>113
h_lim_1 = data[0]<136
l_lim_2 = data[0]>136
h_lim_2 = data[0]<1000
fit_data_1_x = data[0][l_lim_1*h_lim_1]
fit_data_2_x = data[0][l_lim_2*h_lim_2]
res1 = least_squares(residuals, x0=p0_1, args=(fit_data_1_x, data[1][l_lim_1*h_lim_1]), loss='huber')
res2 = least_squares(residuals, x0=p0_2, args=(fit_data_2_x, data[1][l_lim_2*h_lim_2]), loss='huber')

x_range1 = np.linspace(fit_data_1_x[0], fit_data_1_x[-1], 100)
x_range2 = np.linspace(fit_data_2_x[0], fit_data_2_x[-1], 100)
plt.plot(x_range1, fct(*res1.x, *[x_range1]))
plt.plot(x_range2, fct(*res2.x, *[x_range2]))
plt.scatter(data[0], data[1])
plt.show()
