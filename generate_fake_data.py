#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from astropy.table import Table
import astropy.constants as const
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['text.usetex'] = True
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.size'] = 11

def get_oscillation_array(M, R, l_max, n_max):
    oscillation_array = np.zeros((n_max+1, l_max+1))
    delta_nu = M**(1/2) * R**(-3/2) * 134.9 / (M_sun**(1/2) * R_sun**(-3/2))
    for i in range(l_max+1):
        for j in range(n_max+1):
            oscillation_array[j, i] = delta_nu * (j + i/2 + 1.45) - 2.5 * i * (i + 1)
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
T_0 = 436 # K
delta_T = 1250 # K
visibility = np.array([1, 1.5, 0.5, 0.04])
duration = 3.154e7 # s
print('Please input the filename of your table or press ENTER for default:')
filename = input('[Target_Selection.csv] >>> ')
if filename == '':
    filename = 'Target_Selection.csv'
t = Table.read(filename)
print('Would you like to run through all stars in the table or select one at random?')
while True:
    selection = input('[all/id/random/sun] >>> ')
    if selection == 'all':
        IDs = t['ID']
        break
    if selection == 'id':
        print('To add a star to the list, input its ID and press ENTER. When you are finished, press ENTER again:')
        IDs = []
        while True:
            query = input('>>> ')
            if query == '':
                break
            elif np.count_nonzero(t['ID'] == int(query)) == 1:
                IDs.append(int(query))
                print('Added to list')
            else:
                print('Star not found')
        break
    if selection == 'random':
        IDs = [np.random.choice(t['ID'])]
        break
    if selection == 'sun':
        IDs = ['Sun']
        L = L_sun
        M = M_sun
        R = R_sun
        g = g_solar
        T_eff = T_eff_sun
        break
for ID in IDs:
    print('Now generating: {}'.format(ID))
    if selection != 'sun':
        index = np.where(t['ID'] == ID)
        L = float(10**t['logL'][index] * L_sun)
        M = float(t['Mass'][index] * M_sun)
        R = float(np.sqrt((G * (t['Mass'][index] * M_sun)) / ((10**t['logg'][index]) * 0.01)))
        g = float(10**t['logg'][index] * 0.01)
        T_eff = float(10**t['logTe'][index])
    oscillation_array = get_oscillation_array(M, R, 3, 39)
    T_red = 8907 * (L / L_sun)**-0.093 # K
    beta = (1 - np.exp((T_eff - T_red) / delta_T))
    amplitude = 2.1 * beta * L * M_sun / L_sun / M * (T_eff_sun / T_eff)**2 # ppm
    nu_max = nu_max_sun * (g / g_solar) * np.sqrt(T_eff_sun / T_eff)
    fwhm = 0.66 * nu_max**0.88 # μHz
    sigma = fwhm / (2 * np.sqrt(2 * np.log(2)))
    cadence = 25 # s
    if nu_max + 2 * fwhm > 20000:
        print('WARNING! This star has oscillations exceeding the observational frequency limit. Would you like to simulate using the fast camera instead? [Note: This will generate 10 times more data!]')
        while True:
            warning = input('[y/n] >>> ')
            if warning == 'y':
                cadence = 2.5 # s
                break
            if warning == 'n':
                break
    lower_limit = 1e6 / duration # μHz
    upper_limit = 1e6 / (2 * cadence) # μHz
    x = np.linspace(lower_limit, upper_limit, round(upper_limit / lower_limit - 1)) # μHz
    y = nu_max / nu_max_sun
    sigma_c = sigma_c_sun * (L / L_sun) * ((M_sun / M)**1.5) * ((T_eff_sun / T_eff)**2.75) * np.sqrt(y)
    tau_c = tau_c_sun / y # Ms
    P_granulation = 4 * sigma_c**2 * tau_c / (1 + (2 * np.pi * x * tau_c)**2)
    alpha = 2.95 * y + 0.39 
    gamma_alpha = 3.08 * y + 3.32
    Wdip = 4637 * y - 141
    nu_dip = 2984 * y + 60 
    delta_gamma_dip = -0.47 * y + 0.62
    plt.figure(figsize=(7, 4))
    y_combined = np.zeros_like(x)
    for i in range(oscillation_array.shape[1]):
        y = np.zeros_like(x)
        for j in range(oscillation_array.shape[0]):
            if nu_max < 3900:
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
    plt.ylim(0)
    plt.show()
    np.save('{}.npy'.format(ID), np.vstack((x, y_combined + P_granulation)).T)
