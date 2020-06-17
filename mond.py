import numpy as np
import matplotlib.pyplot as plt
from numpy import sin,cos,pi,sqrt,exp
from numpy.linalg import norm
from random import uniform,gauss
from time import time

start = time()
def euler_bound(x,v):
    x_k = x
    for i in range(n_particles):
        x[i] += v[i]*dt
        for j in range(n_particles):
            if(i!=j):
                a = norm(f(x_k[i],x_k[j]))
                muv = mu(a/a_0)
                #print(a,muv)
                v[i] += (f(x_k[i],x_k[j])/muv)*dt
        if (norm(x[i]) > R) or (norm(x[i]) < -R):
            v[i] = -v[i]
            x[i] += v[i]*dt
    return x,v
def euler(x,v):
    x_k = x
    for i in range(n_particles):
        x[i] += v[i]*dt
        for j in range(n_particles):
            if(i!=j):
                a = norm(f(x_k[i],x_k[j]))
                muv = mu(a/a_0)
                #print(a,muv)
                v[i] += (f(x_k[i],x_k[j])/muv)*dt
    return x,v
def symplectic(x,v):
    for i in range(n_particles):
        x[i] += v[i]*dt
        for j in range(n_particles):
            if(i!=j):    
                v[i] += f(x[i],x[j])*dt
    return x,v
def init_two():
    x1 = ([R*cos(omega*t0),R*sin(omega*t0)])
    x2 = -np.copy(x1)
    v1 = ([omega*x1[1],omega*x1[0]])
    v2 = -np.copy(v1)
    x = np.array([x1,x2])
    v = np.array([v1,v2])
    return x,v
def get_uniform_coordinates():
    x = np.zeros((n_particles,d))
    for i in range(n_particles):
        x[i] = ([uniform(-R,R)*cos(omega*t0),uniform(-R,R)*sin(omega*t0)])
    return x
def get_init_velocities():
    v = np.zeros((n_particles,d))
    for i in range(n_particles):
        v[i] = ([omega*x[i,1],omega*x[i,0]])
    return v
def f(xi,xj):
    rij = xj-xi
    return (G*m*(rij))/(norm(rij)+epsilon)**3
def mu(s):
    return s/sqrt(1+s**2)

#parameter
m = 1. #kg
R = 2. #m
G = 6.67 #m/s^2
omega = sqrt((G*m)/(4*R**3)) #velocities
epsilon = 0.001
d = 2 #dimension
n_particles = 10 #particles
t0 = 0.
t = 5*2.0*pi/omega
dt = 0.01
N = np.int(np.floor(t/dt))
scale = 20.0
a_0 =  1e-8
print(t,N)
#initial condition
#x,v = init_two()
x = get_uniform_coordinates()
v = get_init_velocities()
#main loop
#print(x,v)
plt.plot(x[:,0],x[:,1], 'ro')
for k in range(N):
    t = k*dt
    x,v = euler_bound(x,v)
    #plt.plot(xe[:,0],xe[:,1], 'b.')
    #plt.xlim(right=scale,left=-scale)
    #plt.ylim(top=scale,bottom=-scale)
    #plt.axes(aspect='equal')
    plt.plot(x[:,0],x[:,1], 'b.')
#filename='./figures/plot.png'
#plt.savefig(filename)
print("Time for running code :", time()-start, "seconds")
plt.show()