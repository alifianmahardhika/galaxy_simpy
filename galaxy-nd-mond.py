#function

import matplotlib.pyplot as plt
from numpy import sin,cos,pi,sqrt,exp,floor,zeros
from numpy.random import normal,uniform
from numpy.linalg import norm
#from random import uniform
from time import time
#from mpmath import sech
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

start = time()
def euler_mond(x,v):
    x_k = x
    a = zeros((n_particles,d))
    for i in range(n_particles):
        x[i,0] += v[i,0]*dt
        x[i,1] += v[i,1]*dt
        for j in range(n_particles):
            a[j] = [beta(np.abs(f(x_k[j,0]))/a_0)*f(x_k[j,0]),beta(np.abs(f(x_k[j,1]))/a_0)*f(x_k[j,1])]
            v[j,0] += (a[j,0]*dt) 
            v[j,1] += (a[j,1]*dt)


def euler_newton(x,v):
    x_k = x
    av = zeros((n_particles,d))
    for i in range(n_particles):
        av[i] = [f(x_k[i,0],x_k[i,0],x_k[i,1]),f(x_k[i,1],x_k[i,0],x_k[i,1])]
        x[i,0] += v[i,0]*dt
        x[i,1] += v[i,1]*dt
        v[i,0] += (av[i,0]*dt)
        v[i,1] += (av[i,1]*dt)

def symplectic(x,v):
    x_k = x
    for i in range(n_particles):
        v[i,0] += (f(x_k[i,0])*dt)
        v[i,1] += (f(x_k[i,1])*dt)
        x[i,0] += v[i,0]*dt
        x[i,1] += v[i,1]*dt

def radius():
    r = zeros((n_particles,d))
    i = 0
    while(i<n_particles):
        r1 = uniform(-R,R)
        r2 = uniform(-R,R)
        if(abs(r1**2+r2**2)<(R)**2):
            r[i] = ([r1,r2])
            i+=1
        else:
            i=i
    return r

def omega():
    Nn = n_particles+1
    w = zeros((Nn-1,d))
    A = zeros((Nn-1,d))
    r = radius()
    for j in range(1,Nn-1):
        A[j] = [G * (m/(8*np.abs(sin(pi*j/Nn)))+M_center),G * (m/(8*np.abs(cos(pi*j/Nn)))+M_center)]
        w[j] = [(A[j,0]/(np.sqrt(r[j,0]**2+r[j,1]**2)+epsilon)**3), (A[j,1]/(np.sqrt(r[j,0]**2+r[j,1]**2)+epsilon)**3)]
    return w

def clamp(x, lower, upper):
    return max(min(x, upper), lower)


    
def get_init_coordinates2():
    x = zeros((n_particles,d))
    Nn=n_particles+1
    for j in range(Nn-1): 
        t = j*dt 
        x1 = (r[j,0])*cos(w[j,0]*t+(2*pi*j/Nn))
        x2 = (r[j,1])*sin(w[j,1]*t+(2*pi*j/Nn))
        x[j] = ([clamp(x1,-R,R),clamp(x2,-R,R)])
            
    return x

def mu(s):
    return s/sqrt(1+s**2)

def beta(q):
    return sqrt((1+sqrt(1+(4/q**2)))/2)

def get_init_velocities2():
    v = zeros((n_particles,d))
    for i in range(n_particles):
        v[i] = [(-w[i,0])*x[i,1],(w[i,1])*x[i,0]]
    return v

def f(x,xa,xb):
    Nn = n_particles+1
    x_ba = xb-xa
    for j in range(1,n_particles-1):
        A = G*(m*(1/(8*np.abs(sin(pi*j/Nn))))+M_center)
        r_eps = (norm(x_ba)+epsilon)
        f = -A/r_eps**3*x
    return f
    
n_particles = 180 #particles
d = 2 #dimension
m = 10e11/n_particles #[MO]
M_total = m*n_particles
M_center = 18*m
R = 2.9 #[kpc]
G = 13.34*10e-11 #[kpc^3 MO^-1 gy^-2]
epsilon = 0.8
T = 1.5
dt = 0.01
t = 0. #
N = int(floor(T/dt))
scale = 4
a_0 =  2e-2
v_max= 3.41

#initial condition
r  = radius()
w = omega()
x = get_init_coordinates2()
v = get_init_velocities2()

#main loop
for k in range(N):
    euler_newton(x,v)
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    plt.xlim(scale,-scale)
    plt.ylim(scale,-scale)
    ax.scatter(x[:,0], x[:,1], zs=0, zdir='z',c='c')
    no = str(k)
    filename='./figures_new/epsilon1_8/plotnd_clamps'+no+'.png'
    plt.savefig(filename)
    plt.close()
    #np.savetxt('./data/myfile_nd_epsilon1_8'+no+'.csv', np.c_[r,x,v],delimiter=',') #fmt='%.18g',
print("Time for running ", N, "iteration :", time()-start, "seconds")
plt.hist2d(x[:,0], x[:,1], cmap=plt.cm.jet)
plt.colorbar()
filename='./figures_new/epsilon1_8/plotnd_distribution'+no+'.png'
plt.savefig(filename)
plt.close()