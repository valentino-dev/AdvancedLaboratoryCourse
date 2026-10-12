import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from olib import *

plt.rcParams.update({
    "font.size": 12,
})

data_aer = np.loadtxt('data/el_res_air.txt')[2:]
times_aer = data_aer[:, 0]
voltage_aer = data_aer[:, 1]

data_atc = np.loadtxt('data/ambient_thermal_conductivity.txt')[2:]
times_atc = data_atc[:, 0]
voltage_atc = data_atc[:, 1]

data_vtc = np.loadtxt('data/vacuum_thermal_conductivity.txt')[2:]
times_vtc = data_vtc[:, 0]
voltage_vtc = data_vtc[:, 1]

data_ver = np.loadtxt('data/el_res_vac.txt')[2:]
times_ver = data_ver[:, 0]
voltage_ver = data_ver[:, 1]



current_aer = np.arange(-5, 6, 1)*2e-2
current_ver = np.arange(-5, 6, 1)*2e-2
current_vtc = np.array([0, 1, 2, 3, 4, 5]) * 1e-3
current_atc = np.array([0, 1, 2, 3, 4, 5]) * 1e-3

ranges_aer = np.array([[60, 90], [100, 135], [145, 180], [185, 220], [230, 265], [0, 50], [275, 310], [320, 365], [375, 410], [420, 445], [455, 520]])
mask_aer = np.array([(times_aer > rang[0]) * (times_aer < rang[1]) for rang in ranges_aer])

ranges_vtc = [[0, 100], [240, 300], [450, 510], [650, 700], [1150, 1250], [1450, 1550]]
mask_vtc = np.array([(times_vtc > rang[0]) * (times_vtc < rang[1]) for rang in ranges_vtc])

ranges_atc = [[0, 100], [275, 450], [580, 680], [850, 950], [1100, 1200], [1400, 1600]]
mask_atc = np.array([(times_atc > rang[0]) * (times_atc < rang[1]) for rang in ranges_atc])

ranges_ver = np.array([[25, 70], [80, 140], [150, 210], [230, 300], [320, 420], [550, 660], [430, 520], [680, 750], [790, 880], [900, 970], [1000, 1100]])
mask_ver = np.array([(times_ver > rang[0]) * (times_ver < rang[1]) for rang in ranges_ver])



aer = np.array([voltage_aer[mask].mean() for mask in mask_aer])
aer_err = np.array([voltage_aer[mask].std() for mask in mask_aer])

ver = np.array([voltage_ver[mask].mean() for mask in mask_ver])
ver_err = np.array([voltage_ver[mask].std() for mask in mask_ver])



atc = np.array([voltage_atc[mask].mean() for mask in mask_atc])
atc_err = np.array([voltage_atc[mask].std() for mask in mask_atc])

vtc = np.array([voltage_vtc[mask].mean() for mask in mask_vtc])
vtc_err = np.array([voltage_vtc[mask].std() for mask in mask_vtc])

d_m=0.111e-3
d_m_err=0.0015e-3
b_m=15.28e-3
b_m_err=0.09e-3
l_TC_m = 49.8e-3
l_TC_m_err = 0.5e-3
l_R_m = 48.8e-3
l_R_m_err = 2.0e-3
r_PH_m=1.67e3

def model_ver(x, a, b):
    return a*x+b
def model_vtc(x, a, b):
    return a*x**2+b

table_atc = Table(current_atc, np.abs(current_atc*0), atc, atc_err, xlabel=r'$I_\text{PH}$', ylabel=r'$U_\text{TH}$', xunit='A', yunit='V', title='Thermal conductivity at ambient pressure.')
fig, ax = plt.subplots()
ax, model_return = plotData(ax, table_atc, model=model_vtc)
ax.set_title(None)
fig.savefig('images/atc.pdf', dpi=500, bbox_inches='tight')
plt.clf()
table_atc.saveAsPDF('./latex/', 'atc_data.tex')


table_vtc = Table(current_vtc, np.abs(current_vtc*0), vtc, vtc_err, xlabel=r'$I_\text{PH}$', ylabel=r'$U_\text{TH}$', xunit='A', yunit='V', title='Thermal conductivity in vacuum.')
fig, ax = plt.subplots()
ax, model_return = plotData(ax, table_vtc, model=model_vtc)
ax.set_title(None)
fig.savefig('images/vtc.pdf', dpi=500, bbox_inches='tight')
plt.clf()
table_vtc.saveAsPDF('./latex/', 'vtc_data.tex')

table_aer = Table(current_aer, np.abs(current_aer*0), aer, aer_err, xlabel='$I_R$', ylabel='$U_R$', xunit='A', yunit='V', title='Electric resistivity measurement at ambient pressure.')
fig, ax = plt.subplots()
ax, aer_returns = plotData(ax, table_aer, model=model_ver)
ax.set_title(None)
fig.savefig('images/aer.pdf', dpi=500, bbox_inches='tight')
plt.clf()
table_aer.saveAsPDF('./latex/', 'aer_data.tex')

table_ver = Table(current_ver, np.abs(current_ver*0), ver, ver_err, xlabel='$I_R$', ylabel='$U_R$', xunit='A', yunit='V', title='Electric resistivity measurement in vacuum.')
fig, ax = plt.subplots()
ax, ver_returns = plotData(ax, table_ver, model=model_ver)
ax.set_title(None)
fig.savefig('images/ver.pdf', dpi=500, bbox_inches='tight')
plt.clf()
table_ver.saveAsPDF('./latex/', 'ver_data.tex')


mask_ = [True]*11
mask_[5] = False

print('aer2')
def calc_resistivity(aer, aer_err, aer_returns, current_aer, mask):
    aer = ((aer-aer_returns[0][1])/current_aer)[mask_]
    aer_err = (((aer_err/current_aer)**2+(aer_returns[1][1,1]/current_aer)**2)**(1/2))[mask_]
    current_aer = current_aer[mask_]

    aeR = aer*d_m*b_m/l_R_m
    aeR_err = ((aeR/aer*aer_err)**2+(aeR/d_m*d_m_err)**2+(aeR/b_m*b_m_err)**2+(aeR/l_R_m*l_R_m_err)**2)**(1/2)
    return aeR, aeR_err, current_aer

aeR, aeR_err, current_aer = calc_resistivity(aer, aer_err, aer_returns, current_aer, mask_)
veR, veR_err, current_ver = calc_resistivity(ver, ver_err, ver_returns, current_ver, mask_)

F = np.average(aeR, weights=1/aeR_err**2)
F_err = np.sqrt(1/(aeR_err**(-2)).sum())
print('aer', roundToError(F, F_err))

F = np.average(veR, weights=1/veR_err**2)
F_err = np.sqrt(1/(veR_err**(-2)).sum())
print('ver', roundToError(F, F_err))

table_aer = Table(current_aer, current_aer*0, aeR, aeR_err, xlabel='$I_R$', ylabel=r'$\rho$', xunit='A', yunit='$\\Omega$m')
fig, ax = plt.subplots()
ax = plotData(ax, table_aer, model=model_ver)
fig.savefig('images/aer2.pdf', dpi=500, bbox_inches='tight')
plt.clf()

table_ver = Table(current_ver, current_ver*0, veR, veR_err, xlabel='$I_R$', ylabel=r'$\rho$', xunit='A', yunit='$\\Omega$m')
fig, ax = plt.subplots()
ax = plotData(ax, table_ver, model=model_ver)
fig.savefig('images/ver2.pdf', dpi=500, bbox_inches='tight')
plt.clf()

S = 5.9e-5
def calc_thermal_cunductivity(atc, atc_err, current_atc):
    u_ph = current_atc*r_PH_m
    tc = u_ph*current_atc*S*l_TC_m/atc/(d_m*b_m)
    tc_err = ((tc/l_TC_m*l_TC_m_err)**2+(tc/atc*atc_err)**2+(tc/d_m*d_m_err)**2+(tc/b_m*b_m_err)**2)**(1/2)
    return tc, tc_err

atC, atC_err = calc_thermal_cunductivity(atc, atc_err, current_atc)
vtC, vtC_err = calc_thermal_cunductivity(vtc, vtc_err, current_vtc)

print('atc2')
table_atc = Table(current_atc[1:], current_atc[1:]*0, atC[1:], atC_err[1:], xlabel='$I_{{PH}}$', ylabel='$\\kappa$', xunit='A', yunit=r'$\frac{\text{W}}{\text{mK}}$')
fig, ax = plt.subplots()
ax = plotData(ax, table_atc, model=model_ver)
fig.savefig('images/atc2.pdf', dpi=500, bbox_inches='tight')
plt.clf()

plt.close('all')
plt.scatter(current_vtc[1:], current_vtc[1:]**2/vtc[1:])
plt.scatter(current_atc[1:], current_atc[1:]**2/atc[1:])
#plt.scatter(current_vtc[1:], vtc[1:])
#plt.scatter(current_atc[1:], atc[1:])
plt.show()

table_vtc = Table(current_vtc[1:], current_vtc[1:]*0, vtC[1:], vtC_err[1:], xlabel='$I_{{PH}}$', ylabel='$\\kappa$', xunit='A', yunit=r'$\frac{\text{W}}{\text{mK}}$')
fig, ax = plt.subplots()
ax = plotData(ax, table_vtc, model=model_ver)
fig.savefig('images/vtc2.pdf', dpi=500, bbox_inches='tight')
plt.clf()

def calc_temp_grad(atc, atc_err):
    atc = atc - atc[0]
    temp_grad = atc/S/l_TC_m
    temp_grad_err = ((atc_err/S/l_TC_m)**2+(atc/S/l_TC_m**2*l_TC_m_err)**2)**(1/2)
    sqrt_tmp_grad = np.sqrt(temp_grad)
    sqrt_tmp_grad_err = temp_grad**(-3/2)*temp_grad_err
    print(temp_grad)
    print(sqrt_tmp_grad)
    print(sqrt_tmp_grad_err)
    return sqrt_tmp_grad, sqrt_tmp_grad_err,

atg, atg_err = calc_temp_grad(atc, atc_err)
vtg, vtg_err = calc_temp_grad(vtc, vtc_err)
table_atc = Table(current_atc, current_atc*0, atg, atg_err, xlabel='$I_{{PH}}$', ylabel=r'$\sqrt{\nabla T}$', xunit='A', yunit=r'$\sqrt{\frac{K}{m}}$')
'''


fig, ax = plt.subplots()
ax = plotData(ax, table_atc)#, model=model_vtc)
#ax.set_yscale('log')
fig.savefig('images/atg.pdf', dpi=500, bbox_inches='tight')
plt.clf()
'''

table_vtc = Table(current_vtc, current_vtc*0, vtg, vtg_err, xlabel='$I_{{PH}}$', ylabel=r'$\sqrt{\nabla T}$', xunit='A', yunit=r'$\sqrt{\frac{K}{m}}$')
fig, ax = plt.subplots()
ax = plotData(ax, table_atc, label='Ambient pressure', color='red', connect=True)#, model=model_vtc)
ax = plotData(ax, table_vtc, label='Vacuum', color='blue', connect=True)#, model=model_vtc)
#ax.set_yscale('log')
fig.savefig('images/sqrtplot.pdf', dpi=500, bbox_inches='tight')
plt.clf()
