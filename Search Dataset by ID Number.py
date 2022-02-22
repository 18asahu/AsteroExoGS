# -*- coding: utf-8 -*-
"""
Created on Tue Feb 15 14:23:49 2022

@author: clent
"""

import numpy as np
from astropy.table import Table
from astropy.table import QTable
from astropy.table import Column
import astropy
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure

t = Table.read(r'C:\Users\clent\Documents\Group Studies\K_bytempandgrav.csv')

print(t.colnames)

ID = 639894
mask = (t['ID'] < (ID+1)) & (t['ID'] > (ID-1))
print(t['logg'][mask])
print(t['logTe'][mask])