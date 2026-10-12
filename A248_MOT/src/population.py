import numpy as np
import matplotlib.pyplot as plt
import olib

# TODO: Change all values to own
d_lense = 2.5e-2
d_lense_err =  0.1e-2
#p_meas = 445e-9
p_meas = 125e-9
#p_meas_err = 14e-9
p_meas_err = 10e-9
#p_back = 280e-9
p_back = 16e-9
#p_back_err = 10e-9
p_back_err = 2e-9
d = 7e-2
d_err = 0.5e-2
p_tot_cool = 7.96e-3
p_tot_cool_err = 0.17e-3

def prnt(v, err):
    print(*olib.roundToError(v, err))


p_fluo = p_meas - p_back
p_fluo_err = (p_meas_err**2 + p_back_err**2)**(1/2)
print(*olib.roundToError(p_fluo, p_fluo_err))

p_tot = p_fluo*4*3.141*d**2/(3.141*(d_lense/2)**2)
p_tot_err = ((p_tot/p_fluo*p_fluo_err)**2+(p_tot/d*2*d_err)**2+(p_tot/d_lense*2*d_lense_err)**2)**(1/2)
print(*olib.roundToError(p_tot, p_tot_err))

v_cool = 384.230e12
h = 6.626e-34

n_gamma = p_tot/(h*v_cool)
n_gamma_err = p_tot_err/(h*v_cool)
print(n_gamma, n_gamma_err)
prnt(n_gamma, n_gamma_err)



omega = 1
omega_err = 1
i_tot_cool = p_tot_cool/3.141/omega
i_tot_cool_err = ((p_tot_cool_err/3.141/omega)**2+(p_tot_cool/3.141/omega**2*omega_err)**2)**(1/2)

i_sat = 1.66932e-3
i_sat_err = 0.00035e-3
i_ratio = i_tot_cool/i_sat
i_ratio_err = ((i_tot_cool_err/i_sat)**2+(i_tot_cool/i_sat**2*i_sat_err)**2)*(1/2)
gamma = 6.0666e6
gamma_err = 0.0018e6
delta = -7.4
delta_err = 1.6
r = i_ratio*3.141*gamma/(1+i_ratio+4*(delta/gamma)**2)
r_err = (((3.141*gamma/(1 + 4*i_ratio*(delta/gamma)**2)**2) * i_ratio_err)**2 + (((3.141*i_ratio*(1 + 12*i_ratio*delta**2/gamma**2))/(1 + 4*i_ratio*delta**2/gamma**2)**2) * gamma_err)**2 + (((-8*3.141*i_ratio**2*delta/gamma)/(1 + 4*i_ratio*delta**2/gamma**2)**2) * delta_err)**2)**(1/2)
w = 2.68e-3
w_err = 0.15e-3
I = p_tot/3.141/w**2
I_err = ((p_tot_err/3.141/w**2)**2+(p_tot/3.141/w**3*w_err*2)**2)**(1/2)

prnt(I, I_err)
