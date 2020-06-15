import numpy as np
import matplotlib.pyplot as plt
from numpy import sin,cos,pi,sqrt,exp
from numpy.linalg import norm
from random import uniform,gauss
from time import time

start = time()
def f(a,b):
    r = b-a
    r = (r[0]**2+r[1]**2)**(1/2)
    #r = sqrt(np.dot(r,r))
    return (G*m*(b-a))/(r+epsilon)**3

def euler_bound(x,v):
    for i in range(n_particles):
        for j in range(n_particles):
            if(i!=j):    
                x[i] += v[i]*dt
                v[i] += f(x[i],x[j])*dt
                if (norm(x[i]) > R) or (norm(x[i]) < -R):
                    v[i] = -v[i]
                    x[i] += v[i]*dt
    return x,v

def euler(x,v):
    x_k = x
    for i in range(n_particles):
        for j in range(n_particles):
            if(i!=j):    
                x[i] += v[i]*dt
                v[i] += f(x_k[i],x_k[j])*dt
    return x,v

def euler_modif(x,v):
    for i in range(n_particles):
        x[i] += v[i]*dt
        for j in range(n_particles):
            if(i!=j):
                v[i] += f(x[i],x[j])*dt
    return x,v

def symplectic(x,v):
    for i in range(n_particles):
        for j in range(n_particles):
            if(i!=j):    
                v[i] += f(x[i],x[j])*dt
    for i in range(n_particles):
         x[i] += v[i]*dt
    return x,v

def symplectic_bound(x,v):
    for i in range(n_particles):
        for j in range(n_particles):
            if(i!=j):    
                v[i] += f(x[i],x[j])*dt
    for i in range(n_particles):
        x[i] += v[i]*dt
        if (norm(x[i]) > R) or (norm(x[i]) < -R):
            v[i] = -v[i]
            x[i] += v[i]*dt
    return x,v

def get_init_coordinates_modif():
    x=np.zeros((n_particles,d))
    i = 0
    while(i<n_particles):
        x1 = np.random.normal(-R,R)
        x2 = np.random.normal(-R,R)
        if (np.abs(x1**2+x2**2)<R**2):
            x[i] = ([x1,x2])
            i += 1
    return x

def get_uniform_coordinates():
    x = np.zeros((n_particles,d))
    for i in range(n_particles):
        x[i] = ([uniform(-R,R)*cos(omega*t0),uniform(-R,R)*sin(omega*t0)])
    return x

def get_uniform_coordinates_modif():
    x = np.zeros((n_particles,d))
    for i in range(n_particles):
        x[i] = ([uniform(-R,R)*cos(np.random.normal(0,2*pi)),uniform(-R,R)*sin(np.random.normal(0,2*pi))])
    return x

def get_init_velocities():
    v = np.zeros((n_particles,d))
    for i in range(n_particles):
        v[i] = ([-omega*x[i,1],omega*x[i,0]])
    return v

def get_init_velocities_modif():
    v = np.zeros((n_particles,d))
    for i in range(n_particles):
        v[i] = ([(2*omega)*cos(omega),(2*omega)*sin(omega)])
    return v

#parameter
m = 1. #kg
R = 2. #m
G = 6.67 #m/s^2
omega = sqrt((G*m)/(4*R**3)) #velocities
epsilon = 0.001
d = 2 #dimension
n_particles = 30 #particles
t0 = 0.
t = 3.*2.0*pi/omega
dt = 0.01
N = np.int(np.floor(t/dt))
scale = 20.0
print(t,N)

x = get_init_coordinates2()
v = get_init_velocities()
time = str(dt)
for k in range(N):
    t = k*dt
    x,v = symplectic_bound(x,v)
    plt.xlim(right=scale,left=-scale)
    plt.ylim(top=scale,bottom=-scale)
    plt.title('Plot symplectic for dt ='+time)
    plt.xlabel(r'$x_1$', fontsize=12)
    plt.ylabel(r'$x_2$', fontsize=12)
    plt.plot(x[:,0],x[:,1], 'b.')
    no = np.str(k)
    plt.savefig('./folder./symp_bound/bound_30p_symp'+no+'.png')
    plt.close()
#plt.show()