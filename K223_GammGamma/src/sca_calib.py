import numpy as np
import matplotlib.pyplot as plt
import olib
import pandas as pd

print("SCA Calibration")


for mode in ['SCA', 'CFD']:
    for detector in ['C', 'M']:
        #mode = 'SCA' # SCA or CFD
        #detector = 'C' # M(oving) or C(onstant)

        gated = pd.read_csv(f"data/Measuring_2026-03-09_10-02-21_angle_180_{mode}_Gated_{detector}.dat", delimiter="\t")
        ungated = pd.read_csv(f"data/Measuring_2026-03-09_10-02-21_angle_180_{mode}_UnGated_{detector}.dat", delimiter="\t")

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
        plt.savefig(f'images/{mode}_Calib_{detector}.pdf', dpi=500)
