import Alexandrov_Soto as als
import scipy as sc
import numpy as np
import matplotlib.pyplot as plt


lst_sol = []
par = als.arg
t = als.t
y0 = als.y0
f = open('FLC.txt', 'w')
t_begin = 360  # Деленное на 10
t_end = 374  # - 2
j_begin = 5
j_end = 13

for temp in range(t_begin, t_end, 1):
    T = temp / 10
    Q = 3 ** ((T - 20) / 10)
    print(f'T = {T}, Q(T) = {round(Q, 3)}')
    f.write(f'T = {T}, Q(T) = {round(Q, 3)}'+'\n***\n')
    print('***')
    for j in range(j_begin, j_end, 1):
        I_syn = j / 10
        par[0] = I_syn
        par[-1] = 3 ** ((T - 20) / 10)
        print(f'I_syn = {I_syn}: ', end=' ')
        f.write(f'I_syn = {I_syn};')
        sol = sc.integrate.odeint(als.diff_ur, y0, t, args=(par,), tfirst=True)
        c_v = 0
        mx_v = max(sol[:, 0])
        mx_n = max(sol[:, 1])
        c_n = 0
        kl_v = -10000
        kl_n = -10000
        p = 0
        for k, el in enumerate(sol):
            if mx_v - el[0] < 0.1 and k - kl_v > 4:
                kl_v = k
                c_v += 1

            if mx_n - el[1] < 0.0008 and k - kl_n > 4:
                c_n += 1
                kl_n = k

        if c_n == c_v and c_n > 1:
            print('LC 99% c_v = c_n =', str(c_n))
            f.write('LC 99% c_v = c_n = '+ str(c_n) + ' diff = '+str(abs(c_n-c_v))+'\n')
        elif abs(c_n - c_v) < 2:
            print(f'LC 80% c_v = {str(c_v)}, c_n = {str(c_n)}')
            f.write(f'LC 80% c_v = {str(c_v)}, c_n = {str(c_n)}, diff = {str(abs(c_n-c_v))} \n')
        else:
            print(f'NLC 99% c_v = {str(c_v)}, c_n = {str(c_n)}')
            f.write(f'NLC 99% c_v = {str(c_v)}, c_n = {str(c_n)}, diff = {str(abs(c_n-c_v))} \n')
    print('----------------')

