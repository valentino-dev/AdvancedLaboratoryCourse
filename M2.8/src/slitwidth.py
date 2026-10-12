print('Slit Width')
from uncertainties import ufloat
from uncertainties import unumpy as unp
import matplotlib.pyplot as plt
import pandas as pd
import glob
import numpy as np
from olib import *

width = [(i+1)*0.1 for i in range(10)]
lables = [f'{width[i]:0.1f}' for i in range(10)]
width.pop(5)
width = np.array(width)
lables.pop(5)
spectra = [pd.read_csv(glob.glob(f'data/SlitWidth/{lables[i]}/*.csv')[0], delimiter=';').to_numpy().T for  i in range(9)]

faulty_data = pd.read_csv(glob.glob(f'data/SlitWidth/0.6/*.csv')[0], delimiter=';').to_numpy().T
fig, ax = plt.subplots()
xb, yb = faulty_data
yb_err = (yb+1)**(1/2)
mask = (xb<560)*(xb>500)
xb = xb[mask]
yb = yb[mask]
yb_err = yb_err[mask]
table = Table(xb, xb*0, yb, yb_err, None, 'Raman shift', 'Intensity', r'$\frac{1}{\text{cm}}$', '(1)')
ax = plotData(ax, table)
plt.savefig('images/bad_data.pdf', dpi=500)

p0=[9e2, 3000, 527, 6, 0.2, 1.5]
def model(x, c0, A, x0, gamma, alpha, sigma):
    return c0 + A * (alpha*gamma**2/((x-x0)**2+gamma**2)+(1-alpha)*np.exp(-(x-x0)**2/(2*sigma**2)))

#p0=[9e2, 3000, 527, 6]
#def model(x, c0, A, x0, gamma):
    #return c0 + A * (gamma**2/((x-x0)**2+gamma**2))

#p0=[9e2, 3000, 527, 1.5]
#def model(x, c0, A, x0, sigma):
    #return c0 + A * np.exp(-(x-x0)**2/(2*sigma**2))


def olivero_fwhm(f_G, f_L):
    """
    Computes the FWHM of a Voigt/Pseudo-Voigt profile using the 
    Olivero-Longbothum numerical approximation.

    Parameters
    ----------
    f_G : float
        Full Width at Half Maximum of the pure Gaussian component.
    f_L : float
        Full Width at Half Maximum of the pure Lorentzian component.

    Returns
    -------
    float
        The approximated total FWHM of the combined profile.
    """
    term0 = f_G**5
    term1 = 2.69269 * (f_G**4) * f_L
    term2 = 2.42843 * (f_G**3) * (f_L**2)
    term3 = 4.47163 * (f_G**2) * (f_L**3)
    term4 = 0.07842 * f_G * (f_L**4)
    term5 = f_L**5
    
    f_p = (term0 + term1 + term2 + term3 + term4 + term5) ** 0.2
    return f_p

popts = list()
popts_err = list()
for i, spectrum in enumerate(spectra):
    fig, ax = plt.subplots()
    x, y = spectrum
    mask = (x<560)*(x>500)
    x = x[mask]
    y = y[mask]
    table = Table(x, x*0, y, (y+1)**(1/2), None, 'Raman shift', 'Intensity', r'$\frac{1}{\text{cm}}$', '(1)')
    #plt.scatter(x, y)
    X = np.linspace(x.min(), x.max(), 1000)
    #plt.plot(X, model(X, *p0))
    ax, (popt, pcov, redchi_sq, sig) = plotData(ax, table, model=model, p0=p0, show_p0=False)
    popts.append(popt)
    popts_err.append(np.sqrt(np.diag(pcov)))
    #plt.show()
    #plt.close('all')
popts = np.array(popts)
popts_err = np.array(popts_err)

FWHM_l = unp.uarray(np.abs(2*popts[:,3]), 2*popts_err[:,3])

FWHM_g = unp.uarray(np.abs(2*popts[:,3]*(2*np.log(2))**(1/2)), 2*popts_err[:,3]*(2*np.log(2))**(1/2))

FWHM = olivero_fwhm(FWHM_g, FWHM_l)
plt.close('all')
fig, ax = plt.subplots()
table = Table(width*2, width*0, unp.nominal_values(FWHM), unp.std_devs(FWHM), None, 'Slit Width', 'FWHM', 'mm', r'cm$^{-1}$')
plotData(ax, table)
plt.savefig('images/slitwidth.pdf', dpi=500)


plt.close('all')
fig, ax = plt.subplots()
table = Table(width*2, width*0, popts[:,4], popts_err[:,4], None, 'Slit Width', r'$\alpha$', 'mm', r'(1)')
plotData(ax, table)
plt.savefig('images/slitwidth_alpha.pdf', dpi=500)
