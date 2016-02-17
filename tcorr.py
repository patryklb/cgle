# -*- coding: utf-8 -*-
"""
Created on Tue Jan  5 18:54:24 2016

@author: patlip
"""

import os
import numpy as np
import matplotlib.pyplot as plt

array = np.array([0])
N = 0
for i in os.listdir(os.getcwd() + '/sims'):
    if i.startswith("S"):
        array_d = np.loadtxt('sims/' + i + '/' + 'corr.txt', usecols = (0, ))
        array_real, array_imag = np.loadtxt('sims/' + i + '/' + 'corr.txt', unpack=True, usecols = (1, 2))
        array = array + array_real + 1j * array_imag
        N = N + 1
    else:
        continue
    
array = array / N
error = 2. * np.sqrt(np.var(array))

x = array_d
y = np.abs(array)

plt.plot(x, y, 'k', color='#1B2ACC')
plt.fill_between(x, y-error, y+error,
    alpha=0.2, edgecolor='#1B2ACC', facecolor='#089FFF',
    linewidth=2, linestyle='-', antialiased=True)
plt.ylabel('Correlation function')
plt.xlabel('x')
plt.title('Correlation function for N = ' + str(N) + ' simulations')
plt.show()
