##
##
## Modelling 1A: Failing chemical reactions simulation.
## Made by Frank van der Top and Timothy van der Valk
##
##

import math as m
import numpy as np
import matplotlib.pyplot as plt


a = 2.0
b = 4.5


def du_dt(__u, __v, __t, __tf) -> float:
    if __t <= __tf:
        return a - b * __u + __u**2 * __v - __u
    else:
        return a * m.exp(__tf - __t) - b * __u + __u**2 * __v - __u


def dv_dt(__u, __v, __t, __tf) -> float:
    return b * __u - __u**2 * __v


def simulate_chemical_process(__u0, __v0, __tf, __tend, __dt):
    """ 
    Simulate the chemical process described by the differential equations.
    Return a tuple of the lists of concentrations for U and V, and a list
    """ 
    n = int(round(__tend / __dt))

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
        du = du_dt(u, v, t, __tf)
        dv = dv_dt(u, v, t, __tf)

        # Compute new values and store.
        u += __dt * du
        v += __dt * dv
        u_list[i] = u
        v_list[i] = v

    return (u_list, v_list, t_list)


def plot_spiral_graph(__maxn, __tend, __tfstart, __tfend, __dt):
    """
    Plot __maxn graphs with tf between 0 and __tend. Print the highest value for V
    with the corresponding tf.
    """
    plt.figure("Spiral Graph")
    plt.xlabel('U')
    plt.ylabel('V')

    # Keep track of highest final value for V.
    max_tf = 0
    max_v = 0
    for i in range(__maxn):
        print("PLOT " + str(i))
        tf = __tfstart + __tfend * (float(i) / float(__maxn))
        (u_list, v_list, t_list) = simulate_chemical_process(0, 0, tf, __tend, __dt)
        plt.plot(u_list, v_list)
    
        maxv = v_list[-1]
        if maxv > max_v:
            max_v = maxv
            max_tf = tf

    plt.title("{:d} graphs with tf in [ 0, {:.2f} ]. Max V = {:.5f}, TF = {:.5f}".format(__maxn, __tend, max_v, max_tf))


def plot_time_graph(__maxn, __tend, __dt):
    """
    Plot __maxn graphs with tf between 0 and __tend. Print the highest value for V
    with the corresponding tf.
    """
    plt.figure("Time Graph")
    plt.xlabel('Time (seconds)')
    plt.ylabel('Concentration (mol)')

    # Keep track of highest final value for V.
    max_tf = 0
    max_v = 0
    for i in range(__maxn):
        tf = __tend * (float(i) / float(__maxn))
        (u_list, v_list, t_list) = simulate_chemical_process(0, 0, tf, __tend, __dt)
        plt.plot(t_list, v_list)
    
        maxv = max(v_list) 
        if maxv > max_v:
            max_v = maxv
            max_tf = tf

    plt.title("{:d} graphs with tf between [ 0, {:.2f} ]. Max V = {:.2f}, TF = {:.2f}".format(__maxn, __tend, max_v, max_tf))


def list_graph(__tf, __tend, __dt):
    """
    Print a list of the values T | U V for given simulation.
    """
    u_list, v_list, t_list = simulate_chemical_process(0, 0, __tf, __tend, __dt)
    for i in range(len(t_list)):
        print("{:2.2f}|  {:e}  {:e}".format(t_list[i], u_list[i], v_list[i]))


def determine_optimal_dt():
    tf = 1.0
    tcompare = 2.0
    error_margin = 0.0001

    # Perform high accuracy simulation.

    dt = 0.2
    for i in range(100):
        dt = 0.2 / (2 ** i)

        # dt simulation.
        u_high, v_high, t_high = simulate_chemical_process(0, 0, tf, tcompare, dt)
        # 2 * dt simulation.
        u_low , v_low, t_low = simulate_chemical_process(0, 0, tf, tcompare, 2 * dt)
  
        err = abs(v_high[-1] - v_low[-1])
        if err < error_margin:
            print("The optimal dt is {:e} with i = {:d}. This is 0.2 / 2^{:d}".format(dt, i, i))
            break


def plot_vector_field(__t, __tf):
    density = 30
    xli = []
    yli = []
    duli = []
    dvli = []

    for yi in range(density):
        for xi in range(density):
            x = xi / density * 5.0
            y = yi / density * 5.0
            du = 1.0
            dv = 0.0

            du = du_dt(x, y, __t, __tf)
            dv = dv_dt(x, y, __t, __tf)

            # Normalize arrows.
            dist = m.sqrt(du * du + dv * dv) + 0.01
            du /= dist
            dv /= dist

            xli.append(x)
            yli.append(y)
            duli.append(du)
            dvli.append(dv)

    plt.quiver(xli, yli, duli, dvli)


def plot_end_v(__tfstart, __tfend, ):
    dt = 0.01
    tend = 10.0
    tflist = [] 
    vlist = []

    n = int(round((__tfend - __tfstart) / dt))
    for i in range(n):
        print("SIM " + str(i))
        tf = __tfstart + (i / n) * (__tfend - __tfstart)
        u_sim, v_sim, t_sim = simulate_chemical_process(0, 0, tf, tend, dt)
        vend = v_sim[-1]
        tflist.append(tf)
        vlist.append(vend)

    plt.figure("Final V as function of TF")
    plt.xlabel('TF (seconds)')
    plt.ylabel('Concentration V (in mol)')
    plt.plot(tflist, vlist)
    plt.show() 


plot_end_v(2, 4)
