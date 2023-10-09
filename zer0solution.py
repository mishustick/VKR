import scipy as sc
import numpy as np
import matplotlib.pyplot as plt
import Alexandrov_Soto as al


def ur(y, par):
    V, n = y
    m_inf = 1 / (1 + np.exp(-(V + 33.8) / 5.2))
    h_na = 1 / (1 + np.exp((V + 60.5) / 9.9))
    n_inf = 1 / (1 + np.exp(-(V + 35) / 5))
    tau_n = 68 / (np.exp(-(25 + V) / 15) + np.exp((30 + V) / 20))
    C = n_inf + h_na
    I_syn, C_m, V_na, V_k, V_l, g_na, g_k, g_l, h_k, Q = par
    I_na = g_na * m_inf**3 * (C - n)*(V - V_na)
    I_k = g_k * n**4 * h_k * (V - V_k)
    I_l = g_l * (V - V_l)

    return [(I_syn - I_na - I_k - I_l),
            (n_inf - n) * Q]


def findroot(i, t):
    v0 = -60
    n0 = 1 / (1 + np.exp(-(v0 + 35) / 5))
    x0 = np.array((v0, n0))
    arg = al.arg
    arg[0] = i
    arg[-1] = 3**((t - 20) / 10)
    zero_solution = sc.optimize.fsolve(ur, x0, args=arg)
    return zero_solution


#print(findroot(0, 36))
#arg = al.arg
#arg[0] = 0
#arg[-1] = 3**((36 - 20) / 10)
#print(ur(findroot(0, 36), arg))
