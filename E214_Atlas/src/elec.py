import numpy as np
import matplotlib.pyplot as plt
import olib



print("Electron Momentum-Energy Ratio")

momentum = np.array([-7.55, -3.27, -89.35, -50.35, -28.62, -32.92, -34.30, -36.40, -55.46, -69.17]) #GeV
energy = np.array([56.6, 47.7, 86.5, .6, 28.1, 64.3, .6, 33.8, 48.8, 49.2])

ratio = momentum/energy
print(momentum)
print(energy)
print(ratio)

plt.hist(ratio, 200)
plt.xlabel("ratio")
plt.ylabel("count")
plt.savefig("images/ME-ratio.pdf", dpi=500)


