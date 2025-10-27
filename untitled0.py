# -*- coding: utf-8 -*-
"""
Created on Wed Oct 22 19:29:04 2025

@author: zacha
"""

from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

def crn_odes(t, y, kb=2, ku=0.2):
    
    A, B, C, D, AC, BC = y 
    
    dAdt = -kb*A*C + ku*AC
    dBdt = -kb*B*C + ku*BC
    dCdt = -kb*A*C - kb*B*C 
    dACdt = kb*A*C - ku*AC
    dBCdt = kb*B*C - ku*BC
    dDdt = 10*AC + 5*BC
    
    return [dAdt, dBdt, dCdt, dDdt, dACdt, dBCdt]

tspan = (0, 100)

plt.figure()

tspan = (0, 100)

A = 10
B = 0
C = 100
D = 0
AC = 0
BC = 0

y0 = [A, B, C, D, AC, BC]

sol0 = solve_ivp(crn_odes, tspan, y0, method='RK45')

plt.plot(sol0.t, sol0.y[3].T, label=f'A = {A}, B = {B}')

A = 0
B = 10

y0 = [A, B, C, D, AC, BC]

sol1 = solve_ivp(crn_odes, tspan, y0, method='RK45')

plt.plot(sol1.t, sol1.y[3].T, label=f'A = {A}, B = {B}')

A = 5
B = 5

y0 = [A, B, C, D, AC, BC]

sol2 = solve_ivp(crn_odes, tspan, y0, method='RK45')

plt.plot(sol2.t, sol2.y[3].T, label=f'A = {A}, B = {B}')

plt.legend()
plt.show()
