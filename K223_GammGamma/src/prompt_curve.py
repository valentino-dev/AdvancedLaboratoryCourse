import numpy as np
import olib
import matplotlib.pyplot as plt
from scipy.optimize import least_squares
from scipy.special import erf

print("Prompt Curve")

data = np.loadtxt('data/fast_concidens_counts.dat', delimiter=',').T
data[0] = data[0]-32
#t=10s d1=32ns

data_x = data[0]
data_y = data[1]
data_y_err = np.sqrt(data[1]+1)
print(data_y_err)

def fit_fct(t, A, A0, sig, w, t0):
    return A/2*(1+erf((t-(t0-w/2))/sig)*erf(((t0+w/2)-t)/sig))+A0

def residuals(params, x, y):
    return fit_fct(*params, *[x]) - y

p0 = [80, 1, 5, 34, 2]
#res = least_squares(residuals, x0=p0, args=(data_x, data_y), loss='huber')
plt.errorbar(data_x, data_y, data_y_err, ls='none', marker='d', capsize=10)
olib.fit_function(fit_fct, data_x, data_y, data_y_err, p0, True)
plt.xlabel(r"Offset $\Delta t/ns$")
plt.ylabel(r"Counts")

plt.savefig(f'images/prompt_curve.pdf', dpi=500, bbox_inches="tight")

