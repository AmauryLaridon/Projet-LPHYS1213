# -*- coding: utf-8 -*-

from scipy.integrate import solve_ivp
import matplotlib.pyplot as mpl

Pr = 0.01
f_2_0 = 1
th_1_0 = -0.07
        

def couche_lim(eta, f):
    return (f[1],f[2],-3*f[0]*f[2] + 2*(f[1]**2) - f[3], f[4],-3*Pr*f[0]*f[4])
    


if __name__ == "__main__":
    print("Couche limite de vitesses")
    sol = solve_ivp(couche_lim, (0,30), [0,0,f_2_0,1,th_1_0])
    mpl.plot(sol.t,sol.y[1], label = "f'(η)")
    mpl.plot(sol.t,sol.y[3], label = "θ(η)")
    mpl.legend()
    mpl.show()