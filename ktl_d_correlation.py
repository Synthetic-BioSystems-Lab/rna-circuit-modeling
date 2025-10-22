# -*- coding: utf-8 -*-
"""
Created on Tue Oct 21 16:40:04 2025

@author: zacha

from scipy.stats import linregress

slope, intercept, r, p, stderr = linregress(d_lst, ktl_lst)

x_fit = np.linspace(min(d_lst), max(d_lst), 100)
y_fit = slope * x_fit + intercept

"""

import matplotlib.pyplot as plt
import numpy as np

ktl_lst = [0.075, 0.05, 0.015, 0.0075, 0.0025, 0.001]
d_lst = [20, 15, 10, 5, 1, 0]

p = np.poly1d(np.polyfit(d_lst, ktl_lst, deg=1))
x_fit = np.linspace(0, 20, 100)
y_fit = p(x_fit)

# Plotting
plt.figure()
plt.plot(d_lst, ktl_lst, linestyle='', marker='o', color='#277A4C')
plt.plot(x_fit, y_fit, linestyle='--', color='#32A160')
plt.show()