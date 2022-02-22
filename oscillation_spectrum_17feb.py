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
            oscillation_array[j, i] = delta_nu * (j + i/2 + 1.45) - D * i * (i + 1)
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
# datadir = '/home/rlh/Documents/Group Studies/Data/'
t = Table.read('K_bytempandgrav.csv')

ID = np.random.choice(t['ID'])
print('Star ID: {}'.format(ID))
L = 10**t['ID' == ID]['logL'] * L_sun
M = t['ID' == ID]['Mass'] * M_sun
R = np.sqrt((G*(t['ID' == ID]['Mass'] * M_sun))/((10**t['ID' == ID]['logg'])*0.01))
g = 10**t['ID' == ID]['logg']*0.01
T_eff = 10**t['ID' == ID]['logTe']
T_0 = 436 # K
T_red = 8907 * (L / L_sun)**(-0.093) # K
delta_T = 1250 # K

# Temporary code to plot the Sun!
L = L_sun
M = M_sun
R = R_sun
g = g_solar
T_eff = T_eff_sun

D = 1.5 # μHz

oscillation_array = oscillation_array(M, R, 3, 41)
beta = (1 - np.exp((T_eff - T_red) / delta_T))
amplitude = 2.1 * beta * L * M_sun / L_sun / M * (T_eff_sun / T_eff)**2 # ppm
visibility = np.array([1, 1.5, 0.5, 0.04])

nu_max = nu_max_sun*((g/g_solar)*(T_eff/T_eff_sun)**-0.5)
nu_max_sun = 3090 # μHz
fwhm = 0.66 * nu_max**0.88 # μHz
sigma = fwhm / (2 * np.sqrt(2 * np.log(2)))
x = np.linspace(1, 20000, round((20000)/0.0317)) # μHz
sigma_c = sigma_c_sun*(L/L_sun)*((M_sun/M)**1.5)*((T_eff_sun/T_eff)**2.75)*((nu_max/nu_max_sun)**0.5)
tau_c = (nu_max/nu_max_sun)**(-1) * (250/1e6) # Ms
P_granulation = 4 * sigma_c**2 * tau_c / (1 + (2 * np.pi * x * tau_c)**2)

y = (nu_max/nu_max_sun)

alpha = 2.95*y + 0.39 
gamma_alpha = 3.08*y + 3.32
Wdip = 4637*y - 141
nu_dip = 2984*y + 60 
delta_gamma_dip = -0.47*y + 0.62

linewidthlist = []
if nu_max < 4000:
    for i in range(len(nu)):
        linewidth = np.exp((alpha*np.log(nu[i]/nu_max) + np.log(gamma_alpha)) + ((np.log(delta_gamma_dip))/(1+(((2*np.log(nu[i]/nu_dip))/np.log(Wdip/nu_max))**2))))
        linewidthlist.append(linewidth)
        
else:
    for i in range(len(nu)):
        linewidth1 = np.exp((alpha*np.log(nu[i]/nu_max) + np.log(gamma_alpha)))
        linewidthlist.append(linewidth1)
                           

plt.figure(figsize=(8, 5))
y_combined = np.zeros_like(x)
for i in l:
    y = np.zeros_like(x)
    index= 0
    for j in nu[i]:
        y += gaussian(x, sigma, nu_max) * amplitude**2 * 2 / np.pi / linewidth[index] * visibility[i] * lorentzian(x, j, linewidth[index])
        index+=1
    plt.plot(x, y+P_granulation, label=r'$l$ = {}'.format(i))
    y_combined += y    

plt.axvline(x=nu_max, ls='--', color='gray', label=r'$\nu_{\rm max}$')
plt.plot(x, P_granulation, color='k', label='Granulation')
plt.title('Stellar oscillation spectrum')
plt.xlabel(r'Frequency, $\nu$ [$\mu$Hz]')
plt.xscale("log")
plt.yscale("log")
plt.ylabel(r'PSD [ppm$^2 \mu$Hz$^{-1}$]')
plt.legend()
# plt.xlim(2250, np.max(x))
plt.show()
np.savetxt('oscillation_data_{}.csv'.format(ID), np.vstack((x, y_combined)).T, delimiter=',')
