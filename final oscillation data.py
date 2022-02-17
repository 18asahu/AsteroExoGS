# -*- coding: utf-8 -*-
"""
Created on Thu Feb 17 14:29:12 2022

@author: Megan
"""

import numpy as np
import matplotlib.pyplot as plt
#import matplotlib as mpl
#mpl.rcParams['text.usetex'] = True
#mpl.rcParams['font.family'] = 'serif'
#mpl.rcParams['font.size'] = 11

def nu(delta_nu, n, l):
    nu = []
    for i in l:
        nu.append(delta_nu * (n + i/2 + 1.45) - D * i * (i + 1))
    return nu

def lorentzian(x, x_0, width):
    return 1 / (np.pi * width * (1 + ((x - x_0) / width)**2))

def gaussian(x, sigma, mu):
    return np.exp(-(x - mu)**2 / (2 * sigma**2)) #/ (sigma * np.sqrt(2 * np.pi))

pFilePath = 'K dwarfs new.txt'
L_K = []
M_K = []
R_K = []
T_eff_K = []
g_K =[]
G = 6.67e-11
index = []
#same pFilePath as above
with open(pFilePath) as pFile:
    next(pFile)
    for pColumn in pFile:
            data_column = pColumn.split()
            logl = float(data_column[12])
            mass = float(data_column[11])
            logteff = float(data_column[13])
            grav = float(data_column[14])
            ind = float(data_column[0])
            L_K.append(10**logl * 3.846e26)
            M_K.append(mass*1.989e30)
            T_eff_K.append(10**logteff)
            g_K.append((10**grav)*0.01)
            R_K.append(np.sqrt((G*(mass*1.989e30))/((10**grav)*0.01)))
            index.append(ind)
            
print(index[3])
# Define solar parameters
L_sun = 3.846e26 # W
M_sun = 1.989e30 # kg
R_sun = 6.957e8 # m
T_eff_sun = 5777 # K
sigma_c_sun = 60
tau_c_sun = 250/1000000
nu_max_sun = 3090 # microHz
g_solar = 274 # ms-2

# Define stellar parameters
L = L_K[3]
M = M_K[3]
R = R_K[3]
T_eff = T_eff_K[3]
g = g_K[3]

# Define solar parameters
L_sun = 3.846e26 # W
M_sun = 1.989e30 # kg
R_sun = 6.957e8 # m
T_eff_sun = 5777 # K
sigma_c_sun = 60
nu_max_sun = 3090 # microHz

tau_c_sun = 250/1000000


# Define temperatures
T_0 = 436 # K
T_red = 8907 * (L / L_sun)**(-0.093) # K
delta_T = 1250 # K

n = np.arange(40) # Array of n values to be plotted
l = np.arange(4) # Array of l values to be plotted
D = 1.5 # μHz

# old version of linewidth
# gamma_0 = 1.02 # μHz
# gamma = gamma_0 * np.exp((T_eff - 5777) / T_0) # μHz

delta_nu = M**(1/2) * R**(-3/2) * 135 / (M_sun**(1/2) * R_sun**(-3/2)) # μHz
nu = nu(delta_nu, n, l) # μHz
beta = (1 - np.exp((T_eff - T_red) / delta_T))
amplitude = 2.1 * beta * L * M_sun / L_sun / M * (T_eff_sun / T_eff)**2 # ppm
visibility = np.array([1, 1.5, 0.5, 0.04])

nu_max = 3090 # μHz
fwhm = 0.66 * nu_max**0.88 # μHz
sigma = fwhm / (2 * np.sqrt(2 * np.log(2)))
x = np.linspace(1, nu_max+fwhm, np.int((nu_max+fwhm-1)/0.0317)) # μHz
sigma_c = sigma_c_sun*(L/L_sun)*((M_sun/M)**1.5)*((T_eff_sun/T_eff)**2.75)*((nu_max/nu_max_sun)**0.5)
tau_c = 250/1000000
P_granulation = 4 * sigma_c**2 * tau_c / (1 + (2 * np.pi * x * tau_c)**2)
# P_granulation = 2 / np.pi * 560**2 / 2.3 / (1 + (x/2.3)**2)

# #new method of linewidth
y = (nu_max/3090)

#wrong way round in paper
alpha = 2.95*y + 0.39 
gamma_alpha = 3.08*y + 3.32
Wdip = 4637*y - 141
nu_dip = 2984*y + 60 
if nu_max < 4000:
    delta_gamma_dip = -0.47*y + 0.62
else:
    delta_gamma_dip = 0.001

linewidthlist = []
for i in range(len(nu)):
    linewidth = np.exp((alpha*np.log(nu[i]/nu_max) + np.log(gamma_alpha)) + ((np.log(delta_gamma_dip))/(1+(((2*np.log(nu[i]/nu_dip))/np.log(Wdip/nu_max))**2))))
    linewidthlist.append(linewidth)


plt.figure(figsize=(8, 5))


y_combined = np.zeros_like(x)
for i in l:
    y = np.zeros_like(x)
    index = 0
    for j in nu[i]:
        y += gaussian(x, sigma, nu_max) * amplitude**2 * 2 / np.pi / linewidth[index] * visibility[i] * lorentzian(x, j, linewidth[index])
        index += 1
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
plt.savefig('oscillation_spectrum_36.png')

curve, = plt.plot(x, y_combined + P_granulation, lw=0)

xdata = curve.get_xdata()
ydata = curve.get_ydata()

np.savetxt('oscillation_data_36.csv', np.vstack((xdata, ydata)).T, delimiter=',')