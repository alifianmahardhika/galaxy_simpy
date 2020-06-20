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
def get_init_velocities():
    v = zeros((n_particles,d))
    for i in range(n_particles):
        v[i] = ([omega*cos(0),omega*sin(0)])
    return v
def f(xi,xj):
    rij = xj-xi
    return (m*(rij))/(norm(rij)+epsilon)**3
def init_two():
    x1 = ([R*cos(omega*0),R*sin(omega*0)])
    x2 = -copy(x1)
    v1 = ([omega*x1[1],omega*x1[0]])
    v2 = -copy(v1)
    x = array([x1,x2])
    v = array([v1,v2])
    return x,v
#Global parameter
n_particles = 2 #particles
d = 2 #dimension
m = 10e11/n_particles #[MO]
R = 2.9 #[kpc]
G = 13.34*10e-11 #[kpc^3 MO^-1 gy^-2]
omega = sqrt((G*m)/(4*R**3)) #velocities
epsilon = 1e-5
T = 5
dt = 0.001
N = int(floor(T/dt))
scale = 10.0
#initial condition
x,v = init_two()
#x = get_init_coordinates()
#v = get_init_velocities()
print(x)
#main loop
plt.plot(x[:,0],x[:,1], 'ro')
for k in range(N):
    symplectic(x,v)
    #plt.plot(xe[:,0],xe[:,1], 'b.')
    #plt.xlim(right=scale,left=-scale)
    #plt.ylim(top=scale,bottom=-scale)
    #plt.axes(aspect='equal')
    plt.plot(x[:,0],x[:,1], 'b.')
#filename='./figures/plot.png'
#plt.savefig(filename)
print("Time for running ", N, "iteration :", time()-start, "seconds")
print(x)
plt.show()