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
    x_k = x
    for i in range(n_particles):
        for j in range(n_particles):
            if(i!=j):    
                x[i] += v[i]*dt
                v[i] += f(x_k[i],x_k[j])*dt
    return x,v, print(x_k,'and',x)

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
        if (norm(x[i]) > R_gal) or (norm(x[i]) < -R_gal):
            v[i] = -v[i]
            x[i] += v[i]*dt
    return x,v

def get_init_coordinates_modif():
    x=np.zeros((n_particles,d))
    i = 0
    while(i<n_particles):
        x1 = np.random.normal(-R_gal,R_gal)
        x2 = np.random.normal(-R_gal,R_gal)
        if (np.abs(x1**2+x2**2)<R_gal**2):
            x[i] = ([x1,x2])
            i += 1
    return x

def get_init_velocities_modif():
    v = np.zeros((n_particles,d))
    for i in range(n_particles):
        v[i] = ([(5*omega)*cos(omega2),(5*omega)*sin(omega2)])
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
omega2 = np.random.normal(0,np.int(2*np.pi))

#particle motion
v_max = 3.27 * (kpc/gal_year)
t0 = 0. #second
t = (2.0*pi/omega)
dt = 4e5* gal_year
N_iter = np.floor(t/dt)
N = np.int(np.abs(N_iter))


#array of particles
d = 2 #dimension
n_particles = np.int(20)

x = get_init_coordinates_modif()
v = get_init_velocities_modif()

scale = 1.5*R_gal
time = str(dt)

for k in range(N):
    t = k*dt
    x,v = euler_bound(x,v)
    plt.xlim(right=scale,left=-scale)
    plt.ylim(top=scale,bottom=-scale)
    plt.axes(aspect='equal')
    plt.title('Plot Euler vs Exact and dt ='+time)
    plt.xlabel(r'$x_1$', fontsize=12)
    plt.ylabel(r'$x_2$', fontsize=12)
    plt.plot(x[:,0],x[:,1], 'b.')
    no = np.str(k)
    plt.savefig('./folder./euler_gal_bound/eul_gam'+no+'.png')
    plt.close()
    #plt.show()
    #no = str(k)
    #filename='./figures/'+no+'.png'
    #plt.savefig(filename)
    #plt.close()
    #no = str(k)
    #np.savetxt("output1/data"+no+".csv", np.c_[x,v], delimiter=",")
#print("Time for running code :", time()-start, "seconds")
#plt.show()