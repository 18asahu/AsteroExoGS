#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from astropy.table import Table
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
            oscillation_array[j, i] = delta_nu * (j + i/2 + 1.45) - 1.5 * i * (i + 1)
    return oscillation_array

def lorentzian(x, x_0, width):
    return 1 / (1 + ((x - x_0) / width / 2)**2)

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

# Import the data table, you will need to edit this
datadir = '/home/rlh/Documents/Group Studies/Data/'
t = Table.read(datadir+'K_bytempandgrav.csv')

ID = np.random.choice(t['ID'])
print('Star ID: {}'.format(ID))
index = np.where(t['ID'] == ID)
L = float(10**t['logL'][index] * L_sun)
M = float(t['Mass'][index] * M_sun)
R = float(np.sqrt((G * (t['Mass'][index] * M_sun)) / ((10**t['logg'][index]) * 0.01)))
g = float(10**t['logg'][index] * 0.01)
T_eff = float(10**t['logTe'][index])
T_0 = 436 # K
T_red = 8907 * (L / L_sun)**(-0.093) # K
delta_T = 1250 # K

# Uncomment to plot the solar oscillation spectrum
#L = L_sun
#M = M_sun
#R = R_sun
#g = g_solar
#T_eff = T_eff_sun

oscillation_array = oscillation_array(M, R, 3, 41)
beta = (1 - np.exp((T_eff - T_red) / delta_T))
amplitude = 2.1 * beta * L * M_sun / L_sun / M * (T_eff_sun / T_eff)**2 # ppm
visibility = np.array([1, 1.5, 0.5, 0.04])
nu_max = nu_max_sun*((g / g_solar)*(T_eff / T_eff_sun)**-0.5)
nu_max_sun = 3090 # μHz
fwhm = 0.66 * nu_max**0.88 # μHz
sigma = fwhm / (2 * np.sqrt(2 * np.log(2)))
x = np.linspace(1, 20000, round((20000)/0.0317)) # μHz
y = nu_max / nu_max_sun
sigma_c = sigma_c_sun * (L / L_sun) * ((M_sun / M)**1.5) * ((T_eff_sun / T_eff)**2.75) * (y**0.5)
tau_c = tau_c_sun / y # Ms
P_granulation = 4 * sigma_c**2 * tau_c / (1 + (2 * np.pi * x * tau_c)**2)
alpha = 2.95 * y + 0.39 
gamma_alpha = 3.08 * y + 3.32
Wdip = 4637 * y - 141
nu_dip = 2984 * y + 60 
delta_gamma_dip = -0.47 * y + 0.62

plt.figure(figsize=(8, 5))
y_combined = np.zeros_like(x)
for i in range(oscillation_array.shape[1]):
    y = np.zeros_like(x)
    for j in range(oscillation_array.shape[0]):
        if nu_max < 4000:
            linewidth = np.exp((alpha*np.log(oscillation_array[j, i]/nu_max) + np.log(gamma_alpha)) + ((np.log(delta_gamma_dip))/(1+(((2*np.log(oscillation_array[j, i]/nu_dip))/np.log(Wdip/nu_max))**2))))
        else:
            linewidth = np.exp((alpha*np.log(oscillation_array[j, i]/nu_max) + np.log(gamma_alpha)))
        y += gaussian(x, sigma, nu_max) * amplitude**2 * 2 / np.pi / linewidth * visibility[i] * lorentzian(x, oscillation_array[j, i], linewidth)
    plt.plot(x, y + P_granulation, label=r'$l$ = {}'.format(i))
    y_combined += y
plt.axvline(x=nu_max, ls='--', color='gray', label=r'$\nu_{\rm max} =$'+' {} '.format(round(nu_max))+r'$\mu$Hz')
plt.plot(x, P_granulation, color='k', label='Granulation')
plt.title('Stellar oscillation spectrum [ID: {}]'.format(ID))
plt.xlabel(r'Frequency, $\nu$ [$\mu$Hz]')
plt.ylabel(r'PSD [ppm$^2 \mu$Hz$^{-1}$]')
plt.legend()
plt.xlim(np.min(x), nu_max+2000)
plt.show()
#plt.savefig('oscillation_spectrum.pdf') 
np.save('oscillation_data_{}.npy'.format(ID), np.vstack((x, y_combined)).T)
