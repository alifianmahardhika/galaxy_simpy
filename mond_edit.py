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
                if (norm(x[i]) > R_gal) or (norm(x[i]) < -R_gal):
                    v[i] = -v[i]
                    x[i] += v[i]*dt
    return x,v

def euler(x,v):
    for i in range(n_particles):
        for j in range(n_particles):
            if(i!=j):    
                x[i] += v[i]*dt
                v[i] += f(x[i],x[j])*dt
    return x,v

def symplectic(x,v):
    for i in range(n_particles):
        for j in range(n_particles):
            if(i!=j):    
                v[i] += f(x[i],x[j])*dt
                x[i] += v[i]*dt
    return x,v

def init_two():
    x1 = ([R_gal*cos(omega*t0),R_gal*sin(omega*t0)])
    x2 = -np.copy(x1)
    v1 = ([-omega*x1[1],omega*x1[0]])
    v2 = -np.copy(v1)
    x = np.array([x1,x2])
    v = np.array([v1,v2])
    return x,v

def get_init_coordinates():
    x = np.zeros((n_particles,d))
    for i in range(n_particles):
        if (i==0):
            x[i] = ([R_gal*cos(omega*t0),R_gal*sin(omega*t0)])
        elif (i%2):
            x[i] = ([-i*R_gal*cos(omega*t0),-i*R_gal*sin(omega*t0)])
        elif (i%2==0):
            x[i] = ([i*R_gal*cos(omega*t0),i*R_gal*sin(omega*t0)])
    return x

def get_gauss_coordinates():
    sigma_x = R_gal * 10
    sigma_y = R_gal * 0.1
    x = np.zeros((n_particles,d))
    for i in range(n_particles):
        x[i] = ([gauss(mu=0, sigma=sigma_x),gauss(mu=0, sigma=sigma_y)])
    return x

def get_uniform_coordinates():
    x = np.zeros((n_particles,d))
    for i in range(n_particles):
        x[i] = ([uniform(-R_gal,R_gal)*cos(omega*t0),uniform(-R_gal,R_gal)*sin(omega*t0)])
    return x

def get_init_velocities():
    v = np.zeros((n_particles,d))
    for i in range(n_particles):
        v[i] = ([-omega*x[i,1],omega*x[i,0]])
    return v


#unity
mass_sun = 1.989e30   #kg
n_total = 1e9         #total
kpc = 1.0e20            #meters
gal_year = 1.0e15       #seconds

#parameter
m = 0.03 * mass_sun #kg
R_gal = 5. * kpc #m
epsilon = 0.001
G = 6.67e-11 #m/s^2
omega = sqrt((G*m)/(4*R_gal**3)) #velocities

#particle motion
v_max = 3.27 * (kpc/gal_year)
t0 = 0. #second
t = (2.0*pi/omega)
dt = 3e4* gal_year
N_iter = np.floor(t/dt)
N = np.int(np.abs(N_iter))


#array of particles
d = 2 #dimension
n_particles = np.int(20)

x = get_uniform_coordinates()
v = get_init_velocities()

#plt.plot(x[:,0],x[:,1], 'r*')
for k in range(N):
    t = k*dt
    x,v = euler_bound(x,v)
    plt.plot(x[:,0],x[:,1], 'b.')
    no = str(k)
    filename='./figures/'+no+'.png'
    plt.savefig(filename)
plt.close()
    #no = str(k)
    #np.savetxt("output1/data"+no+".csv", np.c_[x,v], delimiter=",")
#print("Time for running code :", time()-start, "seconds")
#plt.show()