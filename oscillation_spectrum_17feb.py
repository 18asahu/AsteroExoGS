# -*- coding: utf-8 -*-
"""
Created on Thu Feb 17 14:29:12 2022

@author: Megan
"""

import numpy as np
from astropy.table import Table
import astropy.units as u
import astropy.constants as const
import matplotlib.pyplot as plt
import matplotlib as mpl
#mpl.rcParams['text.usetex'] = True
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.size'] = 11

def oscillation_array(M, R, l_max, n_max):
    oscillation_array = np.zeros((n_max+1, l_max+1))
    delta_nu = M**(1/2) * R**(-3/2) * 134.9 / (M_sun**(1/2) * R_sun**(-3/2))
    for i in range(l_max+1):
        for j in range(n_max+1):
            oscillation_array[j, i] = delta_nu * (j + i/2 + 1.45) - D * i * (i + 1)
    return oscillation_array

def lorentzian(x, x_0, width):
    return 1 / (np.pi * width * (1 + ((x - x_0) / width)**2))

def gaussian(x, sigma, mu):
    return np.exp(-(x - mu)**2 / (2 * sigma**2))

G = const.G.value
L_sun = const.L_sun.value # W
M_sun = const.M_sun.value # kg
R_sun = const.R_sun.value # m
T_eff_sun = 5777 # K
sigma_c_sun = 60
tau_c_sun = 250e-6 # Ms
nu_max_sun = 3090 # μHz
g_solar = 274 # ms-2

datadir = '/home/rlh/Documents/Group Studies/Data/'
t = Table.read(datadir+'K_bytempandgrav.csv')
L_K = []
M_K = []
R_K = []
T_eff_K = []
g_K =[]
for i in range(len(t)):
    L_K.append(10**t['logL'][i] * L_sun)
    M_K.append(t['Mass'][i] * M_sun)
    T_eff_K.append(10**t['logTe'][i])
    g_K.append((10**t['logg'][i])*0.01)
    R_K.append(np.sqrt((G*(t['Mass'][i] * M_sun))/((10**t['logg'][i])*0.01)))

L = L_K[3]
M = M_K[3]
R = R_K[3]
T_eff = T_eff_K[3]
g = g_K[3]

T_0 = 436 # K
T_red = 8907 * (L / L_sun)**(-0.093) # K
delta_T = 1250 # K

D = 1.5 # μHz

oscillation_array = oscillation_array(M, R, 3, 41)
beta = (1 - np.exp((T_eff - T_red) / delta_T))
amplitude = 2.1 * beta * L * M_sun / L_sun / M * (T_eff_sun / T_eff)**2 # ppm
visibility = np.array([1, 1.5, 0.5, 0.04])

nu_max = 3090 # μHz
fwhm = 0.66 * nu_max**0.88 # μHz
sigma = fwhm / (2 * np.sqrt(2 * np.log(2)))
x = np.linspace(1, nu_max+fwhm, round((nu_max+fwhm-1)/0.0317)) # μHz
sigma_c = sigma_c_sun*(L/L_sun)*((M_sun/M)**1.5)*((T_eff_sun/T_eff)**2.75)*((nu_max/nu_max_sun)**0.5)
tau_c = 250/1000000
P_granulation = 4 * sigma_c**2 * tau_c / (1 + (2 * np.pi * x * tau_c)**2)

y = (nu_max/3090)

alpha = 2.95*y + 0.39 
gamma_alpha = 3.08*y + 3.32
Wdip = 4637*y - 141
nu_dip = 2984*y + 60 
if nu_max < 4000:
    delta_gamma_dip = -0.47*y + 0.62
else:
    delta_gamma_dip = 0.001

plt.figure(figsize=(8, 5))
y_combined = np.zeros_like(x)
for i in range(oscillation_array.shape[1]):
    y = np.zeros_like(x)
    for j in range(oscillation_array.shape[0]):
        linewidth = np.exp((alpha*np.log(oscillation_array[j, i]/nu_max) + np.log(gamma_alpha)) + ((np.log(delta_gamma_dip))/(1+(((2*np.log(oscillation_array[j, i]/nu_dip))/np.log(Wdip/nu_max))**2))))
        y += gaussian(x, sigma, nu_max) * amplitude**2 * 2 / np.pi / linewidth * visibility[i] * lorentzian(x, oscillation_array[j, i], linewidth)
    plt.plot(x, y + P_granulation, label=r'$l$ = {}'.format(i))
    y_combined += y
plt.axvline(x=nu_max, ls='--', color='gray', label=r'$\nu_{\rm max}$')
plt.plot(x, P_granulation, color='k', label='Granulation')
plt.title('Stellar oscillation spectrum')
plt.xlabel(r'Frequency, $\nu$ [$\mu$Hz]')
plt.xscale("log")
plt.yscale("log")
plt.ylabel(r'PSD [ppm$^2 \mu$Hz$^{-1}$]')
plt.legend()
plt.savefig('oscillation_spectrum.pdf')
np.savetxt('oscillation_data_36.csv', np.vstack((x, y_combined)).T, delimiter=',')