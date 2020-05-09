# -*- coding: utf-8 -*-

from scipy.integrate import solve_ivp
import matplotlib.pyplot as mpl
        

def couche_lim(eta, f):
    return (f[1],f[2],-3*f[0]*f[2] + 2*(f[1]**2) - f[3], f[4],-3*Pr*f[0]*f[4])

def optim(Pr, f_2_0,th_1_0):
    n = 7
    sol = solve_ivp(couche_lim, (0,40), [0,0,f_2_0,1,th_1_0])
    while abs(sol.y[1][-1]) > 0.001:
        if n > 50 :
            print("n > 50")
            break
        
        if sol.y[1][-1] > 0.001:
            f_2_0 *= 1 - 1/(2**(n+1))
            th_1_0 *= 1 - 1/(2**n)
        elif sol.y[1][-1] < -0.001:
            f_2_0 *= 1 + 1/(2**(n+1))
            th_1_0 *= 1 + 1/(2**n)
            
        sol = solve_ivp(couche_lim, (0,40), [0,0,f_2_0,1,th_1_0])
        n += 1
    return sol
    


if __name__ == "__main__":
    print("Couche limite de vitesses")
    Pr = 0.01
    f_2_0 = 1.01
    th_1_0 = -0.06
    sol = optim(Pr, f_2_0,th_1_0)
    mpl.plot(sol.t,sol.y[1], label = "f'(η)")
    mpl.plot(sol.t,sol.y[3], label = "θ(η)")
    mpl.legend()
    mpl.show()