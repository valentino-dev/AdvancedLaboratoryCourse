import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import least_squares

print("Angular Correlation")


data = np.loadtxt('data/pseudo/angular_correlation.dat', delimiter='\t').T

def fct(b, B, A22, A44, theta):
    theta = theta/360*2*3.141
    return B*(1+A22*(3*np.cos(theta)**2-1)/2+A44*(35*np.cos(theta)**4-30*np.cos(theta)**2+3))+b
def residuals(params, x, y):
    return fct(*params, *[x]) - y
p0 = [1000, 9000, 0.1020, 0.0091]
res = least_squares(residuals, x0=p0, args=(data[0], data[1]), loss='huber')
plt.plot(data[0], fct(*p0, *[data[0]]))
plt.plot(data[0], fct(*res.x, *[data[0]]))
plt.scatter(data[0], data[1])
plt.show()

