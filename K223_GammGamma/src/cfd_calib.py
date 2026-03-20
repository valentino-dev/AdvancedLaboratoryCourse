import numpy as np
import matplotlib.pyplot as plt
import olib
import pandas as pd

print("CFD Calibration")

gated = pd.read_csv("data/Measuring_2026-03-09_10-02-21_angle_180_CFD_Gated_C.dat", delimiter="\t")
ungated = pd.read_csv("data/Measuring_2026-03-09_10-02-21_angle_180_CFD_UnGated_C.dat", delimiter="\t")

fig, ax1 = plt.subplots()
ax1.scatter(gated['ch'], gated['N'], c='green', s=2, marker='.', label='Gated')
ax1.scatter(ungated['ch'], ungated['N'], c='red', s=2, marker='.', label='Ungated')


sensitivity = gated['N']/ungated['N']
ax2 = ax1.twinx()
ax2.plot(gated['ch'], sensitivity, c='blue', linewidth=1, label=r'sensitivity $\eta$')

ax1.set_xlabel('Channel')
ax1.set_ylabel('Count N')
ax2.set_ylabel(r'sensitivity $\eta$')

fig.legend()
fig.tight_layout()
plt.savefig('images/CFD_Calib.pdf', dpi=500)
