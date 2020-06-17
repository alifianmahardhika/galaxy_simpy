import numpy as np
import matplotlib.pyplot as plt
from numpy import sin,cos,pi,sqrt,exp
from numpy.linalg import norm
from random import uniform,gauss
from time import time

start = time()

def symplectic(x,v):
    for i in range(n_particles):
        for j in range(n_particles):
            if(i!=j):    
                v[i] += f(x[i],x[j])*dt
                x[i] += v[i]*dt
    return x,v

def init_two():
    x1 = ([R*cos(omega*t0),R*sin(omega*t0)])
    x2 = -np.copy(x1)
    v1 = ([omega*x1[1],omega*x1[0]])
    v2 = -np.copy(v1)
    x = np.array([x1,x2])
    v = np.array([v1,v2])
    return x,v

#parameter
m = 1. #kg
R = 2. #m
G = 6.67 #m/s^2
omega = sqrt((G*m)/(4*R**3)) #velocities
epsilon = 0.001
d = 2 #dimension
n_particles = 2 #particles
t0 = 0.
t = 5*2.0*pi/omega
dt = 1
N = np.int(np.floor(t/dt))
scale = 20.0
print(t,N)
#initial condition
x,v = init_two()
#main loop
def f(xi,xj):
    rij = xj-xi
    return (G*m*(rij))/(norm(rij)+epsilon)**3
def mu(s):
    return s/sqrt(1+s**2)
a_0 = 10**(-8)
acc = 10**(-9)/a_0
print(a_0,acc,mu(acc))
plt.plot(x[:,0],x[:,1], 'ro')
for k in range(10):
    x_k = x
    for i in range(n_particles):
        x[i] += v[i]*dt
        for j in range(n_particles):
            if(i!=j):
                a = norm(f(x_k[i],x_k[j]))
                muv = mu(a/a_0)
                #print(a,muv)
                v[i] += (f(x_k[i],x_k[j])/muv)*dt
    plt.plot(x[:,0],x[:,1], 'b.')
print("Time for running code :", time()-start, "seconds")
plt.show()