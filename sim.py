##
##
## Modelling 1A: Failing chemical reactions simulation.
## Made by Frank van der Top and Timothy van der Valk
##
##

import math as m
import numpy as np
import matplotlib.pyplot as plt

def simulate_chemical_process(__a0, __b0, __u0, __v0, __tf, __tend, __dt):
    """ 
    Simulate the chemical process described by the differential equations.
    Return a tuple of the lists of concentrations for U and V, and a list
    """ 
    n = int(round(__tend / __dt))
    a = __a0
    b = __b0

    u_list = np.zeros(n + 1)
    v_list = np.zeros(n + 1)
    t_list = np.zeros(n + 1)
    u_list[0] = __u0
    v_list[0] = __v0
    for i in range(n+1):
        t_list[i] = __dt * i

    for i in range(1, n + 1):
        u = u_list[i - 1]
        v = v_list[i - 1]
        t = t_list[i - 1]

        # Start reducing supply of A if t > tf.
        du = 0
        if t <= __tf:
            du = a - b * u + u**2 * v - u
        else: 
            du = a * m.exp(-(t - __tf)) - b * u + u**2 * v - u
        dv = b * u - u**2 * v

        # Compute new values and store.
        u += __dt * du
        v += __dt * dv
        u_list[i] = u
        v_list[i] = v

    return (u_list, v_list, t_list)

(u_list, v_list, t_list) = simulate_chemical_process(2.0, 4.5, 0, 0, 15.0, 30.0, 0.01)

plt.figure("concentration U and V")
plt.plot(t_list, u_list)
plt.plot(t_list, v_list)
plt.xlabel('time in t')
plt.ylabel("concentration in mol")

plt.figure('the spiral')
plt.plot(u_list, v_list)
#plt.axis([-10,10, 0, 5])
plt.xlabel('u')
plt.ylabel('y')
plt.show()
