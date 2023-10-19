import scipy as sc
import numpy as np
import matplotlib.pyplot as plt


def diff_ur(t, y, par):
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
    return [(I_syn - I_na - I_k - I_l) / C_m,
            (n_inf - n) * Q / tau_n]


t = np.linspace(0, 100, 2000)
I_s = 0.9
T = 37
Cm = 1
V_Na = 52
V_K = -84
V_L = -63
g_Na = 2.3
g_K = 2.4
g_L = 0.03
h_K = 0.7329
Q = 3**((T - 20) / 10)
#Q = 8.4
arg = np.array((I_s, Cm, V_Na, V_K, V_L, g_Na, g_K, g_L, h_K, Q))
y0 = np.array((-59.109828596342865, 0.109366837166741))

if __name__ == '__main__':
    print(Q)

    sol = sc.integrate.odeint(diff_ur, y0, t, args=(arg, ), tfirst=True)

    fig1 = plt.figure(figsize=(10, 10))
    axv_n = fig1.add_subplot(111)
    axv_n.plot(sol[:, 0], sol[:, 1], 'black')
    #axv_n.plot(SL[:, 0], SL[:, 1], 'grey')
    axv_n.set_xlabel('V, мВ', size=20)
    axv_n.set_ylabel('n', size=20)
    plt.grid()
    axv_n.scatter(y0[0], y0[1])

    fig2 = plt.figure(figsize=(10, 10))
    axt_v = fig2.add_subplot(211)
    axt_v.plot(t, sol[:, 0], color='green')
    plt.grid()
    axt_v.set_xlabel('t, мc', size=20)
    axt_v.set_ylabel('V, мВ', size=20)
    #axt_v.set_ylim([16, 17])
    axt_n = fig2.add_subplot(212)
    axt_n.plot(t, sol[:, 1], color='red')
    axt_n.set_xlabel('t, мc', size=20)
    axt_n.set_ylabel('n', size=20)
    #axt_n.set_ylim([0.75, 0.88])

    plt.grid()
    plt.show()
    print(sol[-1][0], sol[-1][1])
