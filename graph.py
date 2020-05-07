# -*- coding: utf-8 -*-

# commentaire TEST
# putain yes ca fonctionne sa maman !!!

def RRK4(function, t_0, f_0, t_f, dt):
    import matplotlib.pyplot as mpl
    T = [t_0]
    F = [f_0]
    while T[-1] < t_f:
        try:
            t = T[-1]
            x = F[-1]
            K_1 = function(t,x)
            K_2 = function(t + (dt/2), [x_i + (dt/2)*k_i for x_i,k_i in zip(x,K_1)])
            K_3 = function(t + (dt/2), [x_i + (dt/2)*k_i for x_i,k_i in zip(x,K_2)])
            K_4 = function(t + dt, [x_i + dt*k_i for x_i,k_i in zip(x,K_3)])
            T.append(t + dt)
            F.append([x_i + (dt/6)*(k1 + 2*k2 + 2*k3 + k4) for k1,k2,k3,k4,x_i in zip(K_1,K_2,K_3,K_4,x)])
        except OverflowError:
            break
    sol = [T,[f[0] for f in F],[f[3] for f in F]]
    mpl.plot(sol[0],sol[1], label = "f(η)")
    mpl.plot(sol[0],sol[2], label = "θ(η)")
    mpl.legend()
    mpl.show()


Pr = 7
f_2_0 = 0.01
th_1_0 = 0.00001
        

def couche_lim(eta, f):
    return (f[1],f[2],-3*f[0]*f[2] + 2*(f[1]**2) - f[3], f[4],-3*Pr*f[0]*f[4])
    


if __name__ == "__main__":
    print("Couche limite de vitesses")
    RRK4(couche_lim, 0, [0,0,f_2_0,1,th_1_0], 2.4, 0.001)
    
