from numpy import *
import scipy as sc
import sympy as sp
import zer0solution as zs
from Alexandrov_Soto import arg
import time


def df1v(V, n):
    global df1_dv
    return eval(str(df1_dv))


def df1n(V, n):
    global df1_dn
    return eval(str(df1_dn))


def df2v(V, n):
    global df2_dv
    return eval(str(df2_dv))


def df2n(V, n):
    global df2_dn
    return eval(str(df2_dn))


par = arg

t_begin = 3450  # Деленное на 100
t_end = 3740  # - 2
j_begin = 90
j_end = 110

f = open('TT_BIF.txt', 'w')
g_centr = open('TT_centr.txt', 'w')

for j in range(j_begin, j_end, 1):
    I_syn = j / 100
    #print(f'I_syn = {I_syn}')
    f.write(f'I_syn = {I_syn}; ')
    #print('***')
    per = True
    for t in range(t_begin, t_end, 1):
        T = t / 100
        Q = 3 ** ((T - 20) / 10)
        #print(f'T = {T}: ', end=' ')

        n, V = sp.symbols('n V')
        v0, n0 = zs.findroot(I_syn, T)
        pr = arg
        pr[0] = I_syn
        pr[-1] = 3**((T - 20) / 10)
        if abs(max(zs.ur([v0, n0], arg)[0], zs.ur([v0, n0], arg)[1])) < 0.0001:  # Проверка корректности корня
            suc = True
        else:
            suc = False
        #print(f'V0 = {round(v0, 3)}, n0 = {round(n0, 4)}, КН = {suc}', end='; ')
        if per:
            f.write(f'V0 = {round(v0, 3)}; n0 = {round(n0, 4)}; \n')
            per = False
        f.write(f'T = {T};')
        f.write(f'КН = {suc},')

        m_inf = sp.expand(1 / (1 + sp.exp(-(V + 33.8) / 5.2)))
        h_na = sp.expand(1 / (1 + sp.exp((V + 60.5) / 9.9)))
        n_inf = sp.expand(1 / (1 + sp.exp(-(V + 35) / 5)))
        tau_n = sp.expand(68 / (sp.exp(-(25 + V) / 15) + sp.exp((30 + V) / 20)))
        C = sp.expand(n_inf + h_na)
        C_m, V_na, V_k, V_l, g_na, g_k, g_l, h_k = par[1:-1]
        I_na = sp.expand(g_na * m_inf ** 3 * (C - n) * (V - V_na))
        I_k = sp.expand(g_k * n ** 4 * h_k * (V - V_k))
        I_l = sp.expand(g_l * (V - V_l))

        f1 = sp.expand(I_syn - I_na - I_k - I_l)
        f2 = sp.expand((n_inf - n) * Q / tau_n)

        df1_dv = sp.diff(f1, V)
        df1_dn = sp.diff(f1, n)
        df2_dv = sp.diff(f2, V)
        df2_dn = sp.diff(f2, n)
        l = sp.symbols('λ')
        har_mat = sp.Matrix([[df1v(v0, n0) - l, df1n(v0, n0)],
                             [df2v(v0, n0), df2n(v0, n0) - l]])
        dt = sp.det(har_mat)
        rts = sp.roots(dt, l)
        for k, key in enumerate(rts):
            if k == 0:
                l1 = complex(key)
            else:
                l2 = complex(key)
            #print(f'λ{k+1} = {round(key, 4)}', end=', ')
            f.write(f'l{k+1} = {round(key, 4)}; ')
        if l1.imag == 0 and l2.imag == 0 and l1.real < 0 and l2.real < 0:
            #print('УУ')
            f.write('УУ;\n')
        elif l1.imag == 0 and l2.imag == 0 and l1.real > 0 and l2.real > 0:
            #print('НУ')
            f.write('НУ;\n')
        elif abs(l1.real) < 0.0003 and abs(l2.real) < 0.0003:
            #print('Центр')
            if abs(round(l1.real, 4)) == 0 and round(l2.real, 4) == 0:
                g_centr.write(f'{I_syn}, {T}, l1 = {round(l1.imag, 4)}*i, '
                              f'l2 = {round(l2.imag, 4)}*i, ЦЕНТР\n')
        elif l1.imag != 0 and l2.imag != 0 and l1.real < 0 and l2.real < 0:
            #print('УФ')
            f.write('УФ;\n')
        elif l1.imag != 0 and l2.imag != 0 and l1.real > 0 and l2.real > 0:
            #print('НФ')
            f.write('НФ;\n')
        else:
            #print('С')
            f.write('С;\n')
    f.write('-------------\n')
    #print('-------------')

f.close()
g_centr.close()
