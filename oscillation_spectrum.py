#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 3 15:13:26 2022

@author: Raoul
"""

import numpy as np
import matplotlib.pyplot as plt
#import matplotlib as mpl
#mpl.rcParams['text.usetex'] = True
#mpl.rcParams['font.family'] = 'serif'
#mpl.rcParams['font.size'] = 11

def nu(delta_nu, n, l):
    nu = []
    for l in l:
        nu.append(delta_nu * (n + l/2 + 1.45) - D * l * (l + 1))
    return nu

def lorentzian(x, x_0, width):
    return 1 / (np.pi * width * (1 + ((x - x_0) / width)**2))

def gaussian(x, sigma, mu):
    return np.exp(-(x - mu)**2 / (2 * sigma**2)) #/ (sigma * np.sqrt(2 * np.pi))

# Define solar parameters
L_sun = 3.846e26 # W
M_sun = 1.989e30 # kg
R_sun = 6.957e8 # m
T_eff_sun = 5777 # K

# Define stellar parameters
L = L_sun
M = M_sun
R = R_sun
T_eff = T_eff_sun

# Define temperatures
T_0 = 436 # K
T_red = 8907 * (L / L_sun)**(-0.093) # K
delta_T = 1250 # K

n = np.arange(40) # Array of n values to be plotted
l = np.arange(4) # Array of l values to be plotted
D = 1.5 # μHz
gamma_0 = 1.02 # μHz
gamma = gamma_0 * np.exp((T_eff - 5777) / T_0) # μHz
delta_nu = M**(1/2) * R**(-3/2) * 135 / (M_sun**(1/2) * R_sun**(-3/2)) # μHz
nu = nu(delta_nu, n, l) # μHz
beta = (1 - np.exp((T_eff - T_red) / delta_T))
amplitude = 2.1 * beta * L * M_sun / L_sun / M * (T_eff_sun / T_eff)**2 # ppm
visibility = np.array([1, 1.5, 0.5, 0.04])

nu_max = 3090 # μHz
fwhm = 0.66 * nu_max**0.88 # μHz
sigma = fwhm / (2 * np.sqrt(2 * np.log(2)))
x = np.linspace(nu_max-2000, nu_max+2000, 10000) # μHz

plt.figure(figsize=(8, 5))

for l in l:
    y = np.zeros_like(x)
    for n in nu[l]:
        y += gaussian(x, sigma, nu_max) * amplitude**2 * 2 / np.pi / gamma * visibility[l] * lorentzian(x, n, gamma)
    plt.plot(x, y, label=r'$l$ = {}'.format(l))

plt.axvline(x=nu_max, ls='--', color='gray', label=r'$\nu_{\rm max}$')
plt.title('Solar oscillation spectrum')
plt.xlabel(r'Frequency, $\nu$ [$\mu$Hz]')
plt.ylabel(r'PSD [ppm$^2 \mu$Hz$^{-1}$]')
plt.xlim(np.min(x), np.max(x))
plt.ylim(0)
plt.legend()
plt.savefig('oscillation_spectrum.pdf')
