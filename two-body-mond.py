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
        v[i] += a_0*phi_inv(norm(sigmaF)/a_0)*(sigmaF/norm(sigmaF))*dt
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
    return (G*m*rij)/(norm(rij)+epsilon)**3
def init_two():
    x1 = ([R*cos(omega*0),R*sin(omega*0)])
    x2 = -copy(x1)
    v1 = ([omega*x1[1],omega*x1[0]])
    v2 = -copy(v1)
    x = array([x1,x2])
    v = array([v1,v2])
    return x,v
def kinetic_energy():
    sigmaN = 0.0
    for i in range(n_particles):
        sigmaN += 0.5*m*norm(v[i])**2
    return sigmaN
def phi_inv(q):
    return sqrt((q**2+sqrt(1+q**2))/2.0)**(-1)
#Global parameter
n_particles = 2 #particles
d = 2 #dimension
m = 10e11/n_particles #[MO]
R = 2.9 #[kpc]
G = 13.34*10e-11 #[kpc^3 MO^-1 gy^-2]
omega = sqrt((G*m)/(4*R**3)) #velocities
epsilon = 1e-3
T = 100
dt = 0.001
N = int(floor(T/dt))
scale = 30.0
a_0 = 10e-1
#initial condition
x,v = init_two()
#x = get_init_coordinates()
#v = get_init_velocities()
print(x)
#main loop
plt.plot(x[:,0],x[:,1], 'ro')
for k in range(N):
    euler(x,v)
    #print(kinetic_energy())
    #plt.plot(xe[:,0],xe[:,1], 'b.')
    #plt.xlim(right=scale,left=-scale)
    #plt.ylim(top=scale,bottom=-scale)
    #plt.axes(aspect='equal')
    if(k%100==0):
        plt.plot(x[:,0],x[:,1], 'b.')
#filename='./figures/plot.png'
#plt.savefig(filename)
print("Time for running ", N, "iteration :", time()-start, "seconds")
print(x)
plt.show()