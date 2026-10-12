import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from olib import *
import pandas as pd


d_c = 0.1111e-3
d_c_err = 0.0015e-3

b_c = 2.16e-3
b_c_err = 0.05e-3

l_tc = 16.2e-3
l_tc = 0.5e-3

l_r = 12.2e-3
l_r_err = 0.5e-3

r_ph = 1e3
r_ph_err = 0

magnetic_field = np.array([[0, 0], [-2, 2], [-4, 4], [-6, 6]])

ranges_plus = [[0, 26], [26, 30], [30, 34], [34, 38]]
ranges_minus = [[50, 84], [46, 50], [42, 46], [38, 42]]

data_rho = pd.read_csv('data/RHO.mtw', delimiter='\t')
resistance_minus = np.array(data_rho['R-voltage (-)']/data_rho['R-current (-)'])
resistance_plus = np.array(data_rho['R-voltage (+)']/data_rho['R-current (+)'])

area = d_c*b_c
area_err = ((d_c_err*b_c)**2+(d_c*b_c_err)**2)**(1/2)

fig, ax = plt.subplots()
if True:
    for k, (ranges) in enumerate(zip(ranges_plus, ranges_minus)):
        arr = list()
        arr_err = list()
        arr_2 = list()

        for j, rng in enumerate(ranges):
            #for i, resistance_ in enumerate([resistance_minus, resistance_plus]):
            if True:
                temperatures = np.array(data_rho['Temperature'][rng[0]:rng[1]])
                resistances = (resistance_minus[rng[0]:rng[1]]+resistance_plus[rng[0]:rng[1]])/2
                resistances_err = resistances*0

                resistivity = resistances*area/l_r
                resistivity_err = ((resistances_err*area/l_r)**2+(resistances*area_err/l_r)**2+(resistances*area/l_r**2*l_r_err)**2)**(1/2)
                if k == 0:
                    arr_2.append(temperatures)
                    arr.append(resistivity)
                    arr_err.append(resistivity_err)
                else:
                    table = Table(temperatures, temperatures*0, resistivity, resistivity_err, None, '$T$', '$\\rho$', 'K', '$\\Omega$m')

                    if j==0:
                        color='blue'
                    else:
                        color='red'

                    ax = plotData(ax, table, label=f'$B={magnetic_field[k, j]}\\;$T', color=color, fmt=['+', 'x', 'o'][k-1], connect=True)
                    #ax.scatter(temperatures, resistivity*1e12)
        if k==0:
            temperatures = np.concatenate(arr_2)
            resistivity = np.concatenate(arr)
            resistivity_err = np.concatenate(arr_err)
            table = Table(temperatures, temperatures*0, resistivity, resistivity_err, None, '$T$', '$\\rho$', 'K', '$\\Omega$m')
            ax = plotData(ax, table, label=f'$B={magnetic_field[k, 0]}\\;$T')
            fig.savefig(f'images/resistivity{k}.pdf', dpi=500, bbox_inches='tight')
            fig, ax = plt.subplots()
            


fig.savefig(f'images/resistivity_all.pdf', dpi=500, bbox_inches='tight')
fig.clf()




print("Electrical Resistivity at HT")
resistance_plus_ht = resistance_plus[50:84]
resistance_minus_ht = resistance_minus[50:84]

resistances = (resistance_plus_ht+resistance_minus_ht)/2
temperatures = np.array(data_rho['Temperature'][50:84])
mask = temperatures > 100
temperatures = temperatures[mask]
resistances = resistances[mask]

resistances_err = resistances*0
resistivity = resistances*area/l_r
resistivity_err = ((resistances_err*area/l_r)**2+(resistances*area_err/l_r)**2+(resistances*area/l_r**2*l_r_err)**2)**(1/2)
table = Table(temperatures, temperatures*0, resistivity, resistivity_err, None, '$T$', '$\\rho$', 'K', '$\\Omega$m')
fig, ax = plt.subplots()
def model(x, a, b):
    return x*a+b
plotData(ax, table, model=model)
fig.savefig(f'images/resistivity_ht.pdf', dpi=500, bbox_inches='tight')
plt.clf()

resistance_minus_average, resistance_plus_average = [np.array([res[rng] for rng in ranges_plus[1:][0]]).mean() for res in [resistance_minus,resistance_plus]]




resas = [[[[resis for resis in res[rng[0]:rng[1]]] for rng in [ranges_plus[1:], ranges_minus[1:]][i]] for i in range(2)] for res in [resistance_plus, resistance_minus]]
resas = np.array(resas)
resas4 = resas
resas = resas.mean(axis=(0, 1, 2))

resas2 = np.array([[resis for rng in [ranges_plus[0], ranges_minus[0]] for resis in res[rng[0]:rng[1]]] for res in [resistance_plus, resistance_minus]])
resas = resas2.mean(axis=0)
#resas = np.concatenate((resas, resas2))


#resas = resas - rrr
temperatures = np.concatenate((np.array(data_rho['Temperature'])[ranges_plus[0][0]:ranges_plus[0][1]], np.array(data_rho['Temperature'])[ranges_minus[0][0]:ranges_minus[0][1]]))

mask = temperatures < 41
temperatures = temperatures[mask]
resas = resas[mask]
resistances = resas
resistances_err = resas*0
resistivity = resistances*area/l_r
resistivity_err = ((resistances_err*area/l_r)**2+(resistances*area_err/l_r)**2+(resistances*area/l_r**2*l_r_err)**2)**(1/2)
def model2(x, a, b):
    return a+b*x**5
table = Table(temperatures, temperatures*0, resistivity, resistivity_err, None, '$T$', '$\\rho$', 'K', '$\\Omega$m')
fig, ax = plt.subplots()
print('lt')
_, [popt, pcov, _] = plotData(ax, table, model=model2, p0=[1, 2])
fig.savefig(f'images/resistivity_lt.pdf', dpi=500, bbox_inches='tight')
plt.clf()



r0 = popt[0]
r0_err = np.sqrt(pcov[0, 0])
x, x_err = temperatures, temperatures*0
y, y_err = resistivity-r0, (resistivity_err**2+r0_err**2)**(1/2)
print(roundToError(r0, r0_err))
resistances = (resistance_plus[-1]+resistance_minus[-1])/2
t1 = (resistance_plus[0]+resistance_minus[0])/2
print(resistances/t1)
print(resistances, r0)
resistances_err = resistances*0
resistivity = resistances*area/l_r
resistivity_err = ((resistances_err*area/l_r)**2+(resistances*area_err/l_r)**2+(resistances*area/l_r**2*l_r_err)**2)**(1/2)
rrr = resistivity/r0
rrr_err = ((resistivity_err/r0)**2+(resistivity/r0**2*r0_err)**2)**(1/2)
print(rrr_err)
print(roundToError(rrr, rrr_err))
table = Table(x, x_err, y, y_err, None, '$T$', '$\\rho$', 'K', '$\\Omega$m')
fig, ax = plt.subplots()
def model3(x, a):
    return a*x**5
plotData(ax, table, model=model3)
#plt.errorbar(x, y, y_err)
ax.set_xscale('log')
ax.set_yscale('log')
#plt.show()
print('rrr plot')
fig.savefig(f'images/resistivity_rrr.pdf', dpi=500, bbox_inches='tight')
plt.clf()


B0ranges = np.array([0, 10, 20, 25])
bfield = np.array([-2, -4, -6, 2, 4, 6, 0])
resas4 = resas4.reshape((2, 6, 4))
#print(resas4[0,3])
#exit()
resas4 = np.concatenate((resas4, np.array([resistance_plus, resistance_minus])[:, B0ranges].reshape(2, 1, 4)), axis=1)
temp = np.array([10, 20, 30, 40])

fig, ax = plt.subplots()
for i in range(resas4.shape[2]):
    resistances = resas4[:,:,i]
    resistances = resistances.mean(0)
    print(resistances.shape)
    resistances_err = resistances * 0
    resistivity = resistances*area/l_r
    resistivity_err = ((resistances_err*area/l_r)**2+(resistances*area_err/l_r)**2+(resistances*area/l_r**2*l_r_err)**2)**(1/2)
    #for k in range(resistivity.shape[0]):
    if True:
        #table = Table(bfield, bfield*0, resistivity[k], resistivity_err[k], None, '$B$', r'$\rho$', 'T', '$\\Omega$m')
        print(bfield.shape)
        print(resistivity.shape)
        args = bfield.argsort()
        table = Table(bfield[args], bfield*0, resistivity[args], resistivity_err[args], None, '$B$', r'$\rho$', 'T', '$\\Omega$m')
        '''
        if k == 1:
            color = 'red'
        else:
            color = 'blue'
        '''
        #plotData(ax, table, label=f'$T={temp[i]}\\;$K ({'$+$' if k==1 else '$-$'})', color=color)
        ax = plotData(ax, table, label=f'$T={10*(i+1):.0f}\;$K', color=['#D81B60', '#1E88E5', '#FFC107', '#004D40'][i], connect=True)
    #ax.grid()

fig.savefig(f'images/resistivity_T.pdf', dpi=500, bbox_inches='tight')
plt.close('all')
resas4 = resas4.mean(0).T
resas4 = np.concatenate([resas4[:,:-1].reshape(4,2,3).mean(1),resas4[:,-1].reshape(4,1)], axis=1)
bfield = np.array([2, 4, 6, 0])
resistances = resas4
resistances_err = resas4*0
resistivity = resistances*area/l_r
resistivity_err = ((resistances_err*area/l_r)**2+(resistances*area_err/l_r)**2+(resistances*area/l_r**2*l_r_err)**2)**(1/2)


plt.clf()
fig.clf()
plt.close('all')
y = resistivity/resistivity[:, -1].reshape(4, 1)-1
y_err = ((resistivity_err/resistivity[:, -1].reshape(4, 1))**2+(resistivity/resistivity[:, -1].reshape(4, 1)**2*resistivity_err[:, -1].reshape(4, 1))**2)**(1/2)
x = bfield.reshape(1, 4)/resistivity[:, -1].reshape(4, 1)
fig, ax = plt.subplots()
colors = ['green', 'blue', 'red', 'yellow']
tables = [Table(x[i], x[i]*0, y[i], y_err[i], None, r'$\frac{B}{\rho_0}$', r'$\frac{\Delta\rho}{\rho_0}$', r'$\frac{T}{\Omega\text{m}}$', '(1)') for i in range(4)]
x_scaling = min([table.x_scaling for table in tables])
y_scaling = max([table.y_scaling for table in tables])
ax.grid()

def model4(x, a):
    return x*a

F = list()
F_err = list()
for i in range(4):
    print(i)
    tables[i].change_scale(x_scaling, y_scaling)
    _, [popt, pcov, _] = plotData(ax, tables[i], label=f'$T={temp[i]}\\;$K', color=colors[i], model=model4, show_band=False, fit_color=colors[i])
    F.append(popt[0])
    F_err.append(np.sqrt(pcov[0,0]))

fig.savefig(f'images/Kohler.pdf', dpi=500, bbox_inches='tight')
F = np.array(F)
F_err = np.array(F_err)
print(F_err)
F = np.average(F, weights=1/F_err**2)
#F_err = np.average(F_err, weights=1/F_err**2)
F_err = np.sqrt(1/(F_err**(-2)).sum())
print(F_err)
print(roundToError(F, F_err))


