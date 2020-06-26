import matplotlib.pyplot as plt
from numpy import sin,cos,pi,sqrt,exp,floor,zeros,copy,array
from numpy.random import normal
from numpy.linalg import norm
from random import uniform
from time import time

start = time()
def euler(x,v):
    for i in range(n_particles):
        sigmaF = zeros(2)
        for j in range(n_particles):
                if(i!=j): 
                    sigmaF += f(x[i],x[j])
        x[i] += v[i]*dt
        v[i] += G*sigmaF*dt
def symplectic(x,v):
    for i in range(n_particles):
        sigmaF = zeros(2)
        for j in range(n_particles):
                if(i!=j): 
                    sigmaF += f(x[i],x[j])
        v[i] += G*sigmaF*dt
        x[i] += v[i]*dt
def f(xi,xj):
    rij = xj-xi
    return (m*rij)/(norm(rij)+epsilon)**3
def init_two():
    x1 = ([R*cos(omega*0),R*sin(omega*0)])
    x2 = -copy(x1)
    v1 = ([omega*x1[1],omega*x1[0]])
    v2 = -copy(v1)
    x = array([x1,x2])
    v = array([v1,v2])
    return x,v
def kinetic_energy():
    ke = 0.0
    for i in range(n_particles):
        ke += 0.5*m*norm(v[i])**2
    return ke
def potential_energy():
    pe = 0.0
    for i in range(n_particles):
        for j in range(n_particles):
            if(1<=i<j<=n_particles):
                rij = x[j]-x[i]
                pe += (G*m*m)/(norm(rij))**3
    return pe
#Global parameter
n_particles = 2 #particles
d = 2 #dimension
m = 10e11/n_particles #[MO]
R = 2.9 #[kpc]
G = 13.34*10e-11 #[kpc^3 MO^-1 gy^-2]
omega = sqrt((G*m)/(4*R**3)) #velocities
#omega = normal(0,2*pi) #velocities
epsilon = 1e-3
T = 10
dt = 0.001
N = int(floor(T/dt))
scale = 30.0
#initial condition
x,v = init_two()
#x = get_init_coordinates()
#v = get_init_velocities()
print(x)
#main loop
#plt.plot(x[:,0],x[:,1], 'ro')
en_total = zeros(N)
k_array = zeros(N)
for k in range(N):
    symplectic(x,v)
    en_total[k] = kinetic_energy()-potential_energy()
    k_array[k] = k*dt
    #plt.plot(xe[:,0],xe[:,1], 'b.')
    #plt.xlim(right=scale,left=-scale)
    #plt.ylim(top=scale,bottom=-scale)
    #plt.axes(aspect='equal')
    #if(k%10==0):
        #plt.plot(x[:,0],x[:,1], 'b.')
plt.xlabel("T [kpc]")
plt.ylabel("Total Energy")
plt.title("Total Energy of Symplectic Newton Dynamics")
plt.plot(k_array,en_total, 'g-')
#filename='./figures/plot.png'
#plt.savefig(filename)
print("Time for running ", N, "iteration :", time()-start, "seconds")
print(x)
plt.show()