# -*- coding: utf-8 -*-
"""
Created on Wed Mar  2 11:11:58 2022

@author: FvTop
"""
import math as m
import numpy as np
import matplotlib.pyplot as plt

dt = 0.01
end = 30.0
tf=10

def simulate_supplied():
    n = int(round(end / dt))
    a = 2.0
    b = 4.5

    u_list = np.zeros(n + 1)
    v_list = np.zeros(n + 1)
    u_list[0] = 0
    v_list[0] = 0
    
    t_list = np.zeros(n + 1)
    for i in range(n+1):
        t_list[i] = dt * i

    for i in range(1, n + 1):
        u = u_list[i - 1]
        v = v_list[i - 1]
        t = t_list[i - 1]
        du=0
        if t<=tf or True:
            du = a - b * u + u**2 * v - u
        else: 
            du = 2*m.exp(-t) - b*u + u**2 * v - u
        
        dv = b * u - u**2 * v
        u += dt * du
        v += dt * dv

        u_list[i] = u
        v_list[i] = v

    return (u_list, v_list, t_list)

u_list, v_list, t_list = simulate_supplied()

plt.figure("concentration U and V")
plt.plot(t_list, u_list)
plt.plot(t_list, v_list)
plt.xlabel('time in t')
plt.ylabel("concentratoin in mol")

plt.figure('the spiral')
plt.plot(u_list, v_list)
#plt.axis([-10,10, 0, 5])
plt.xlabel('u')
plt.ylabel('y')
plt.show()