import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import least_squares

print("Prompt Curve")

data = np.loadtxt('data/pseudo/prompt_curve.dat', delimiter='\t').T

def gaussian_fct(sig, mu, x):
    return 1/(sig*(2*3.141)**(1/2))*np.exp(-(x-mu)**2/(2*sig**2))

def smered_saw(sig, mu, width, hight, x):
    l_limit, h_limit = int(mu-width/2), int(mu+width/2)
    i = np.arange(l_limit, h_limit)
    if type(x)==np.ndarray:
        y = np.array([(gaussian_fct(sig, 0, i-xs) * hight).sum() for xs in x])
    else:
        y = (gaussian_fct(sig, 0, i-x) * hight).sum()
    return y

data_x = np.arange(data.shape[1])
data_y = data[1]

fit_x = data_x
fit_y = np.array([smered_saw(3, 22, 10, 1100, x) for x in fit_x])

def fit_fct(sig, mu, width, hight, y0, x):
    return smered_saw(sig, mu, width, hight, x) + y0

def residuals(params, x, y):
    return fit_fct(*params, *[x]) - y

p0 = [3, 24, 17, 1100, 100]
res = least_squares(residuals, x0=p0, args=(data_x, data_y), loss='huber')

plt.plot(data_x, fit_fct(*res.x, *[data_x]))
plt.scatter(data_x, data_y)
plt.show()

