import matplotlib.pyplot as plt
from numpy import sin,cos,pi,sqrt,exp,floor,zeros
from numpy.random import normal
from numpy.linalg import norm
from random import uniform
from time import time

start = time()
def euler_bound(x,v):
    for i in range(n_particles):
        x[i] += v[i]*dt
        for j in range(n_particles):
            if(i!=j):
                a = norm(f(x[i],x[j]))
                muv = mu(a/a_0)
                #print(a,muv)
                v[i] += (f(x[i],x[j])/muv)*dt
        if (norm(x[i]) > R) or (norm(x[i]) < -R):
            v[i] = -v[i]
            x[i] += v[i]*dt
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
def symplectic(x,v):
    for i in range(n_particles):
        x[i] += v[i]*dt
        for j in range(n_particles):
            if(i!=j):    
                v[i] += f(x[i],x[j])*dt
def get_init_coordinates():
    x = zeros((n_particles,d))
    i = 0
    while(i<n_particles):
        x1 = normal(-R,R)
        x2 = normal(-R,R)
        if(abs(x1**2+x2**2)<R**2):
            x[i] = ([x1,x2])
            i+=1
        else:
            i=i
    return x
def get_distribution():
    x = get_init_coordinates()
    rho = zeros((n_particles,d))
    i = 0
    while(i<n_particles):
        #rho[i] = [(M_total/(pi*R**2))*exp(x[i,0]/R)*sech((x[i,1]/R))**2,(M_total/(pi*R**2))*exp(x[i,0]/R)*sech((x[i,1]/R))**2]
        i += 1
    return rho,x

def get_init_velocities():
    v = zeros((n_particles,d))
    for i in range(n_particles):
        v[i] = ([-(1/sqrt(rho[i,0]))*omega*sin(omega*t),(1/sqrt(rho[i,1]))*omega*cos(omega*t)])
    return v
def f(xi,xj):
    rij = xj-xi
    return (G*m*(rij))/(norm(rij)+epsilon)**3
def mu(s):
    return s/sqrt(1+s**2)


#Global parameter
n_particles = 10 #particles
d = 2 #dimension
m = 10e11/n_particles #[MO]
M_total = m*n_particles
R = 2.9 #[kpc]
G = 13.34*10e-11 #[kpc^3 MO^-1 gy^-2]
omega = normal(0,2*pi) #velocities
epsilon = 1e-8
T = 5
dt = 0.01
t = 0. #step
N = int(floor(T/dt))
scale = 7.0
a_0 =  1e-10
#initial condition
#x,v = init_two()
rho,x = get_distribution()
plt.plot(x[:,0],x[:,1], 'r.')
v = get_init_velocities()
#main loop
for k in range(N):
    t = k*dt
    euler(x,v)
    #plt.plot(xe[:,0],xe[:,1], 'b.')
    plt.xlim(right=scale,left=-scale)
    plt.ylim(top=scale,bottom=-scale)
    #plt.axes(aspect='equal')
    plt.plot(x[:,0],x[:,1], 'b.')
#filename='./figures/plot.png'
#plt.savefig(filename)
print("Time for running ", N, "iteration :", time()-start, "seconds")
plt.show()