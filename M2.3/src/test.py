from olib import *
import numpy as np
import matplotlib.pyplot as plt

def model(x, a ,b):
    return a*x**2 + b

x = np.arange(0, 10)*10.
x += np.random.rand(len(x))
table = Table(x, x*0, x**2, 10*x)

fig, ax = plt.subplots()
ax = plotData(ax, table, model=model)
ax.legend()
plt.show()
