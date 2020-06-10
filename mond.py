import numpy as np
import matplotlib.pyplot as plt
from numpy import sin,cos,pi,sqrt,exp
from random import uniform,gauss

def f(a,b):
    r = b-a
    r = (r[0]**2+r[1]**2)**(1/2)
    #r = sqrt(np.dot(r,r))
    return (G*m*(b-a))/(r+epsilon)**3

def euler_bound(x,v):
    x_k = x
    for i in range(n_particles):
        for j in range(n_particles):
            if(i!=j):    
                x[i] += v[i]*dt
                v[i] += f(x_k[i],x_k[j])*dt
            if (abs(x[i,0]) > scale) or (abs(x[i,0] < -scale)):
                v[i] = -v[i]
                x[i] += v[i]*dt
            if (abs(x[i,1]) > scale) or (abs(x[i,1] < -scale)):
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

def get_init_coordinates():
    x = np.zeros((n_particles,d))
    for i in range(n_particles):
        if (i==0):
            x[i] = ([R*cos(omega*t0),R*sin(omega*t0)])
        elif (i%2):
            x[i] = ([-i*R*cos(omega*t0),-i*R*sin(omega*t0)])
        elif (i%2==0):
            x[i] = ([i*R*cos(omega*t0),i*R*sin(omega*t0)])
    return x

def get_gauss_coordinates():
    sigma_x = R * 10
    sigma_y = R * 0.1
    x = np.zeros((n_particles,d))
    for i in range(n_particles):
        x[i] = ([gauss(mu=0, sigma=sigma_x),gauss(mu=0, sigma=sigma_y)])
    return x

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

#parameter
m = 1. #kg
R = 2. #m
G = 6.67 #m/s^2
omega = sqrt((G*m)/(4*R**3)) #velocities
epsilon = 0.1
d = 2 #dimension
n_particles = 10 #particles
t0 = 0.
t = 5*2.0*pi/omega
dt = 100
N = np.int(np.floor(t/dt))
scale = 20.0
print(t,N)
#initial condition
#x,v = init_two()
x = get_uniform_coordinates()
v = get_init_velocities()
#main loop
#print(x,v)
plt.plot(x[:,0],x[:,1], 'ro')
for k in range(10):
    t = k*dt
    x,v = euler(x,v)
    #plt.plot(xe[:,0],xe[:,1], 'b.')
    #plt.xlim(right=scale,left=-scale)
    #plt.ylim(top=scale,bottom=-scale)
    #filename='./figures/fig'+str(k)+'.png'
    #plt.savefig(filename)
    #plt.axes(aspect='equal')
    plt.plot(x[:,0],x[:,1], 'b.')
plt.show()