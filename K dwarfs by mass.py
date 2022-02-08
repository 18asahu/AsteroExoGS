# -*- coding: utf-8 -*-
"""
Created on Thu Feb  3 14:51:58 2022

@author: clent
"""
import numpy as np
from astropy.table import Table
from astropy.table import Column
import astropy
import matplotlib.pyplot as plt


t = Table.read(r'C:\Users\clent\Documents\Group Studies\platosim_SPF.csv', delimiter=' ')
#converted .dat to .csv manually through notebook, delimeter is space
print(t['Mass'])

mass_mask_upper = t['Mass'] < 0.85 #upper mass limit
inter = t[mass_mask_upper] #removes high mass stars
mass_mask_lower = inter['Mass'] > 0.65 #lower mass limit
final = inter[mass_mask_lower] #removes low mass stars
print(final)
print(np.max(final['Mass']))
print(np.min(final['Mass'])) #checking!

final.write(r'C:\Users\clent\Documents\Group Studies\K_bymass.csv')