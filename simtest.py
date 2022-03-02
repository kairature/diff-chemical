import numpy as np
import matplotlib.pyplot as plt

dt = 0.1
end = 20.0

def simulate_supplied():
    n = int(round(end / dt))
    a = 2.0
    b = 4.5

    u_list = np.zeros(n + 1)
    v_list = np.zeros(n + 1)
    t_list = np.zeros(n + 1)
    for i in range(n):
        t_list[i] = dt * i

    for i in range(1, n + 1):
        u = u_list[i - 1]
        v = v_list[i - 1]
        t = t_list[i - 1]

        du = a - b * u + u**2 * v - u
        dv = b * u - u**2 * v
        u += dt * du
        v += dt * dv

        u_list[i] = u
        v_list[i] = v

    return (u_list, v_list, t_list)

u_list, v_list, t_list = simulate_supplied()

plt.plot(t_list, u_list, t_list, v_list)
plt.show()
