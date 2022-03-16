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


def spiral_graph(__maxn, __tend, __dt):
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
        tf = __tend * (float(i) / float(__maxn))
        (u_list, v_list, t_list) = simulate_chemical_process(2.0, 4.5, 0, 0, tf, __tend, __dt)
        plt.plot(u_list, v_list)
    
        maxv = max(v_list) 
        if maxv > max_v:
            max_v = maxv
            max_tf = tf

    plt.title("{:d} graphs with tf in [ 0, {:.2f} ]. Max V = {:.2f}, TF = {:.2f}".format(__maxn, __tend, max_v, max_tf))
    plt.show()


def time_graph(__maxn, __tend, __dt):
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
        (u_list, v_list, t_list) = simulate_chemical_process(2.0, 4.5, 0, 0, tf, __tend, __dt)
        plt.plot(t_list, v_list)
    
        maxv = max(v_list) 
        if maxv > max_v:
            max_v = maxv
            max_tf = tf

    plt.title("{:d} graphs with tf between [ 0, {:.2f} ]. Max V = {:.2f}, TF = {:.2f}".format(__maxn, __tend, max_v, max_tf))
    plt.show()


def list_graph(__tf, __tend, __dt):
    u_list, v_list, t_list = simulate_chemical_process(2.0, 4.5, 0, 0, __tf, __tend, __dt)
    for i in range(len(t_list)):
        print("{:2.2f}|  {:e}  {:e}".format(t_list[i], u_list[i], v_list[i]))

spiral_graph(100, 30.0, 0.1)
