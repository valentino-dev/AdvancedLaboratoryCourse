print('Diffraction Peaks')

import numpy as np
from uncertainties import unumpy as unp
import matplotlib.pyplot as plt
import pandas as pd

#spectrum = pd.read_csv('data/as-1.txt', delimiter=(10, 17))[56:].to_numpy()
b_spectrum = np.genfromtxt('data/bs-1.txt', delimiter=(10, 17))[57:].T
a_spectrum = np.genfromtxt('data/as-1.txt', delimiter=(10, 17))[57:].T
plt.plot(*b_spectrum, c='black')
CuOpeaks = [32.86, 35.78, 38.91, 48.98]
CuOhight = np.array([1000, 1000, 1000, 1000])*0+3e2
CuOlabel= ['(002)', '(111)', '(111)', '']

cuo_2theta_list = [32.5, 35.5, 38.8, 48.7, 58.3, 61.5, 66.2, 68.0, 75.2]
cuo_labels_latex = ['CuO ' + string for string in [
    "$(002)$",
    "$(\\bar{1}11)$",
    "$(111)$",
    "$(\\bar{2}02)$",
    "$(202)$",
    "$(\\bar{1}13)$",
    "$(022)$",
    "$(\\bar{2}20)$",
    "$(\\bar{2}22)$"
]]

# Standard 2-Theta peak positions for Cubic Y2O3
y2o3_2theta_list = [20.5, 29.1, 33.8, 35.9, 39.8, 43.5, 48.5, 53.3, 57.6, 60.4, 71.4, 72.3, 78.5]

# Standard Miller index labels using LaTeX math formatting for Matplotlib
y2o3_labels_latex = ['Y$_3$O$_3$ '+ string for string in [
    "$(211)$",
    "$(222)$",
    "$(400)$",
    "$(\\bar{4}11)$",
    "$(422)$",
    "$(431)$",
    "$(440)$",
    "$(611)$",
    "$(622)$",
    "$(\\bar{4}44)$",
    "$(\\bar{8}00)$",
    "$(822)$",
    "$(\\bar{6}62)$"
]]

# Standard 2-Theta peak positions for Orthorhombic BaCO3 (Witherite)
baco3_2theta_list = [23.9, 24.3, 27.7, 34.1, 34.6, 42.0, 44.1, 44.7, 46.8, 50.4, 52.1]

# Standard Miller index labels using LaTeX math formatting for Matplotlib
baco3_labels_latex = [ 'BaCO$_3$ ' + string for string in [
    "$(111)$",
    "$(020)$",
    "$(002)$",
    "$(022)$",
    "$(112)$",
    "$(130)$",
    "$(200)$",
    "$(041)$",
    "$(220)$",
    "$(202)$",
    "$(132)$"
]]

theta_list = np.concatenate((np.array(y2o3_2theta_list), np.array(cuo_2theta_list), np.array(baco3_2theta_list))) + 0.207
lables_list = np.concatenate((np.array(y2o3_labels_latex), np.array(cuo_labels_latex), np.array(baco3_labels_latex)))
args = theta_list.argsort()
theta_list = theta_list[args]
lables_list = lables_list[args]

#pop_list = []
#pop_list = [0, 5, 10, 14, 15, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32]
pop_list = [0, 5, 10, 14, 15, 19, 20, 21, 22, 24, 25, 26, 27, 29, 30, 31, 32]

for i in pop_list[::-1]:
    theta_list = np.delete(theta_list, i)
    lables_list = np.delete(lables_list, i)

def get_color(string):
    color = ''
    match string[0]:
        case 'C':
            color = 'red'
        case 'Y':
            color = 'green'
        case 'B':
            color = 'blue'
    return color

hight = np.array([656, 390, 170, 350, 320, 290, 407, 260, 270, 187, 400, 250, 200, 245,  200, 200])

#hight = np.zeros(len(theta_list))+200
for i in range(0, len(lables_list)):
    plt.annotate(
        f"{lables_list[i]}", 
        xy=(theta_list[i], hight[i]),          # Coordinate of the peak (arrow tip)
        xytext=(theta_list[i], max(3.0e2+(len(lables_list)-i)*40, hight[i]+20)), # Coordinate of the text label
        #xytext=(theta_list[i], 0), # Coordinate of the text label
        arrowprops=dict(
            color=get_color(lables_list[i]), 
            shrink=0, 
            width=0.2, 
            headwidth=1,
            shrinkA=0,
            shrinkB=0,
        ),
        #horizontalalignment='left',
        #verticalalignment='bottom'
    )

plt.plot([-100], [-100], c='red', label='CuO')
plt.plot([-100], [-100], c='green', label='Y$_3$O$_3$')
plt.plot([-100], [-100], c='blue', label='BaCO$_3$')


plt.ylim((0, 1000))
plt.xticks(np.arange(0, b_spectrum[0].max(), 5))
plt.ylabel('Intensity / (1)')
plt.xlabel('$2\\Theta$ / Degrees')
#plt.grid()
plt.xlim((10, 80))
plt.legend(facecolor="white", framealpha=1, fancybox=False)
plt.savefig('images/before_sintering_spectrum.pdf', dpi=500)


plt.close('all')
plt.clf()



# Standard 2-Theta peak positions for Orthorhombic YBa2Cu3O7 (YBCO)
ybco_2theta_list = [7.7, 15.3, 22.9, 27.8, 32.6, 32.9, 38.6, 40.4, 46.7, 47.5, 58.3, 58.8, 68.2, 68.9]

# Standard Miller index labels using LaTeX math formatting for Matplotlib
ybco_labels_latex = ['' + string for string in [
    "$(001)$",
    "$(002)$",
    "$(003)$",
    "$(102)$",
    "$(013)$",
    "$(110)$",
    "$(104)$",
    "$(113)$",
    "$(020)$",
    "$(200)$",
    "$(123)$",
    "$(213)$",
    "$(026)$",
    "$(220)$"
]]

# High-angle secondary phase peak positions for Cubic Y2O3 (75° to 80°)
y2o3_high_angle_2theta = [77.8, 78.8]

# Corresponding Miller index labels using LaTeX math formatting
y2o3_high_angle_labels = [
    "$(840)$",
    "$(\\bar{6}62)$"
]

theta_list = np.concatenate((np.array(ybco_2theta_list), np.array(y2o3_high_angle_2theta)))
lables_list = np.concatenate((np.array(ybco_labels_latex), np.array(y2o3_high_angle_labels)))
args = theta_list.argsort()
theta_list = theta_list[args]
lables_list = lables_list[args]

hight = [290, 370, 300, 1200, 2150, 400, 560, 560, 400, 690, 420, 270, 450, 300, 300]

pop_list = [0]
for i in pop_list[::-1]:
    theta_list = np.delete(theta_list, i)
    lables_list = np.delete(lables_list, i)



def gab(i):
    if i < 5:
        return 0 
    else:
        return 0

def get_color2(i):
    if i > len(lables_list) - 3:
        return 'green'
    else:
        return 'red'

plt.plot(*(a_spectrum), c='black')
for i in range(0, len(lables_list)):
    plt.annotate(
        f"{lables_list[i]}", 
        xy=(theta_list[i], hight[i]),          # Coordinate of the peak (arrow tip)
        xytext=(theta_list[i] + gab(i), max(6.0e2+(len(lables_list)-i)*100, hight[i]+100)), # Coordinate of the text label
        #xytext=(theta_list[i], 0), # Coordinate of the text label
        arrowprops=dict(
            color=get_color2(i), 
            shrink=0, 
            width=0.2, 
            headwidth=1,
            shrinkA=0,
            shrinkB=0,
        ),
        horizontalalignment='left',
        #verticalalignment='bottom'
    )

plt.xlim((10, 90))
plt.ylim((0, 2400))
plt.ylabel('Intensity / (1)')
plt.xlabel('$2\\Theta$ / Degrees')
plt.plot([-100], [-100], c='red', label='YBa$_2$Cu$_3$O$_7$')
plt.plot([-100], [-100], c='green', label='Y$_3$O$_3$')
plt.legend(facecolor="white", framealpha=1, fancybox=False)
plt.xticks(np.arange(10, 95, 5))
plt.savefig('images/after_sintering_spectrum.pdf', dpi=500)
