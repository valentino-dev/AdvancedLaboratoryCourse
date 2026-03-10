import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import least_squares
from scipy.special import erf

print("Threshold Scan")


data = np.loadtxt('data/pseudo/threshold_scan.dat', delimiter=',').T

#plt.scatter(data[0, 1:], data[1, 1:])
#x = np.linspace(-10, 10, 1000)
#def fct(x, x0, omega, sig, A, A0):
    #return A/2*(1+erf((x-(x0-omega))/sig)*erf(((x0+omega)- x)/sig))+A0

#plt.plot(x, fct(x, 1, 1, 10**-1, 1, 1))

plt.show()
