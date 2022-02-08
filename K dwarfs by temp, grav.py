# -*- coding: utf-8 -*-
"""
Created on Tue Feb  8 15:28:09 2022

@author: clent
"""

import numpy as np
from astropy.table import Table
from astropy.table import Column
import astropy
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure


t = Table.read(r'C:\Users\clent\Documents\Group Studies\platosim_SPF.csv', delimiter=' ')
#converted .dat to .csv manually through notebook, delimeter is space

temp_mask_upper = t['logTe'] < 3.723
inter = t[temp_mask_upper]
temp_mask_lower = inter['logTe'] > 3.594
final = inter[temp_mask_lower]

print(final.colnames)

plt.figure(figsize=(8,8))
plt.scatter(final['logTe'], final['logg'], marker='.')
plt.gca().invert_xaxis()
plt.gca().invert_yaxis()
plt.show()

grav_mask_upper = final['logg'] < 6
mid = final[grav_mask_upper] 
grav_mask_lower = mid['logg'] > 3.99
finalg = mid[grav_mask_lower]

print(np.max(finalg['logg']))
print(np.min(finalg['logg']))

plt.figure(figsize=(8,8))
plt.scatter(finalg['logTe'], finalg['logg'], marker='.')
plt.gca().invert_xaxis()
plt.gca().invert_yaxis()
plt.show()