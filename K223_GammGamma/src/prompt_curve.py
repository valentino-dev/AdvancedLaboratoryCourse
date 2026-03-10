import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import least_squares
from scipy.special import erf

print("Prompt Curve")

data = np.loadtxt('data/fast_concidens_counts.dat', delimiter=',').T
data[0] = data[0]-32
#t=10s d1=32ns

data_x = data[0]
data_y = data[1]

def fit_fct(A, A0, sig, w, t0, t):
    return A/2*(1+erf((t-(t0-w/2))/sig)*erf(((t0+w/2)-t)/sig))+A0

def residuals(params, x, y):
    return fit_fct(*params, *[x]) - y

p0 = [80, 1, 5, 34, 2]
res = least_squares(residuals, x0=p0, args=(data_x, data_y), loss='huber')

plt.plot(data_x, fit_fct(*res.x, *[data_x]))
plt.scatter(data_x, data_y)
plt.show()

