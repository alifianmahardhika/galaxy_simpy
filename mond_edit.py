import numpy as np
import matplotlib.pyplot as plt
from numpy import sin,cos,pi,sqrt,exp
from numpy.linalg import norm
from random import uniform,gauss
from time import time


start = time()

#functions
def f(a,b):
    r = b-a
    r = (r[0]**2+r[1]**2)**(1/2)
    #r = sqrt(np.dot(r,r))
    return (G*m*(b-a))/(r+epsilon)**3

def get_uniform_coordinates():
    x = np.zeros((n_particles,d))
    for i in range(n_particles):
        x[i] = ([uniform(-R_gal,R_gal)*cos(omega*t0),uniform(-R_gal,R_gal)*sin(omega*t0)])
    return x

def get_uniform_coordinates3():
    x = np.zeros((n_particles,d))
    i = 0
    while(i<n_particles):
        x1 = uniform(-R_gal,R_gal)
        x2 = uniform(-R_gal,R_gal)
        if (np.abs(x1**2+x2**2)<R_gal**2):
            x[i] = ([x1,x2])
            i += 1
    return x

def get_init_velocities():
    v = np.zeros((n_particles,d))
    for i in range(n_particles):
        v[i] = ([omega*x[i,1],omega*x[i,0]])
    return v

def get_init_velocities3():
    v = np.zeros((n_particles,d))
    for i in range(n_particles):
        v[i] = ([v_max*cos(omega2),v_max*sin(omega2)])
    return v

def euler_bound2(x,v):
    for i in range(n_particles):
        x[i] += v[i]*dt
        for j in range(n_particles):
            if(i!=j):    
                v[i] += f(x[i],x[j])*dt
                if (norm(x[i]) > R_gal) or (norm(x[i]) < -R_gal):
                    v[i] = -v[i]
                    x[i] += v[i]*dt
    return x,v

def euler_bound(x,v):
    x_k = x
    for i in range(n_particles):
        for j in range(n_particles):
            if(i!=j):    
                x[i] += v[i]*dt
                v[i] += f(x_k[i],x_k[j])*dt
                if (norm(x[i]) > R_gal) or (norm(x[i]) < -R_gal):
                    v[i] = -v[i]
                    x[i] += v[i]*dt
    return x,v

def euler(x,v):
    for i in range(n_particles):
        x[i] += v[i]*dt
        for j in range(n_particles):
            if(i!=j):    
                v[i] += f(x[i],x[j])*dt
    return x,v

#unity
mass_sun = 1.989e30   #kg
n_total = 1e9         #total
kpc = 1.0e20            #meters
gal_year = 1.0e15       #seconds

#parameter
m = 0.03 * mass_sun #kg
R_gal = 5. * kpc #m
G = 6.67e-11 #m/s^2

#particle motion
v_max = 3.27 * (kpc/gal_year)
t0 = 0. #second
t = 7. * gal_year #second
dt = 0.01 * gal_year
#rotate = 2.*np.pi/t
N_iter = np.floor(2*t/dt)
N = np.int(N_iter)
omega = sqrt((G*m)/(4*R_gal**3)) #velocities
omega2 = np.random.normal(0,np.int(2*np.pi))
epsilon = 0.01

#array of particles
d = 2 #dimension
n_particles = np.int(35)
#initial condition
#x,v = init_two()
x = get_uniform_coordinates3()
v = get_init_velocities3()
#omega = sqrt((G*m)/(4*(R_gal**3))
#main loop
#print(x,v)
plt.plot(x[:,0],x[:,1], 'r*')
for k in range(N):
    t = k*dt
    x,v = euler_bound2(x,v)
    #plt.plot(xe[:,0],xe[:,1], 'b.')
    #plt.xlim(right=scale,left=-scale)
    #plt.ylim(top=scale,bottom=-scale)
    #plt.axes(aspect='equal')
    plt.plot(x[:,0],x[:,1], 'b.')
#filename='./figures/plot.png'
#plt.savefig(filename)
    #no = str(k)
    #np.savetxt("output1/data"+no+".csv", np.c_[x,v], delimiter=",")
#print(x[:,1])
print("Time for running code :", time()-start, "seconds")
plt.show()