# -*- coding: utf-8 -*-
"""
Created on Tue May  5 18:45:25 2020

@author: Marina
"""
from numpy import *

def RK4(Xstart,Ustart,Xend,h,f):
  imax = int((Xend-Xstart)/h)
  X = Xstart + arange(imax+1)*h
  U = zeros((size(Ustart), imax+1))
  U[:, 0] = asarray(Ustart)
  for i in range(1,imax+1):
        K1 = asarray(f(U[:, i-1]))
        K2 = asarray(f(U[:, i-1] + (h/2) * K1))
        K3 = asarray(f(U[:, i-1] + (h/2)* K2))
        K4 = asarray(f(U[:, i-1] + h * K3))
        
        U[:, i] = U[:, i-1] + (h/6) * (K1 + 2 * K2 + 2 * K3 + K4)
        
  U = transpose(U)
  return X, U

# We define:
# f = u1; f' = v1; f''= w1
# \u03B8 = u2; \u03B8'= v2
# u = [u1, v1, w1, u2, v2]

# In that case the system to solve is:
# u1' = v1
# v1' = w1
# w1' = -3*u1*w1 + 2*(v1)**2 - u2
# u2' = v2
# v2' = -3*Pr*u1*v2

def f(u):
  du1dt =   u[1]
  dv1dt =   u[2] 
  dw1dt = - 3*u[0]*u[2] + 2*u[1]**2 - u[3]
  du2dt = u[4]
  dv2dt = - 3* 0.01 * u[0]*u[4]
  return array([du1dt,dv1dt,dw1dt,du2dt,dv2dt])

# Shoot method

def shoot(a,b,h,integrator):
  X,U = integrator(0,[0,0,a,1,b],5,h,f)
  return array([U[-1,2], U[-1,4]])

# We solve the problem in the same way we solved the Blasius equation,
# aplaying the bisection method.
# We are supposing that we want to get the an error equal to the tolerance in both cases,
# to calculate alfa1 and alfa2, and that we can choose different intervals.

def blasius(delta1,delta2,nmax,tol,h,integrator):
  delta0 = array([delta1, delta2])
  a = zeros(2)
  b = zeros(2)
  a[0] = delta0[0,0]
  a[1] = delta0[1,0]
  b[0] = delta0[0,1]
  b[1] = delta0[1,1]
  x = zeros(2)
  delta = zeros(2)
  
  print('a1,a2',a[0],a[1])
  fa = shoot(a[0],a[1],h,integrator)
  print('fa', fa)

  fb = shoot(b[0],b[1],h,integrator)
  print('fb', fb)
  n = 1
  delta[0] = (b[0]-a[0])/2 
  delta[1] = (b[1]-a[1])/2
  for i in range (0,2):
   if (fa[i]*fb[i] > 0) :
    return 0, 0,'Bad initial interval :-( for i = %d' % (i) 

  while (abs(delta[0]) >= tol and abs(delta[1]) >= tol and n < nmax) :
      n = n + 1
      for i in range (0,2):
       delta[i] = (b[i]-a[i])/2
       x[i] = a[i] + delta[i]
       print('x%d = %13.7e'%(i,x[i]))
   
   
      fx = shoot(x[0],x[1],h,integrator)
      print('fx', fx)
      print(" x1 = %14.7e (Estimated error %13.7e at iteration %d)" % (x[0],abs(delta[0]),n))
      print(" x2 = %14.7e (Estimated error %13.7e at iteration %d)" % (x[1],abs(delta[1]),n))
      
      for i in range (0,2):
       if (fx[i]*fa[i] > 0) :
        a[i] = x[i]
        fa[i] = fx[i]
       else:
        b[i] = x[i]
        fb[i] = fx[i]
 
  if (n == nmax) :
    return x[0], x[1],'Increase nmax : more iterations are needed :-(' 
  return x[0], x[1],'Convergence observed :-)'

# These are the initial value we have choosen

Xend = 5
h = 0.9
Pr = 0.01
delta1 = [1,3]
delta2 = [2,5]
nmax = 10
tol = 1e-2

a, b, message = blasius(delta1,delta2,nmax,tol,h,RK4)
X,U = RK4(0,[0,0,a,1,b],Xend,h,f)

print('For RK4:') 
print(" === Requested f''(0) = %.4f === %s" % (a,message))
print(" === Requested \u03B8'(0) = %.4f === %s" % (b,message))
print(" === Obtained final value for f'(0) = %13.7e " % U[-1,2])
print(" === Obtained final value for \u03B8(0) = %13.7e " % U[-1,4])


