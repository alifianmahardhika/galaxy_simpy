import matplotlib.pyplot as plt
import numpy as np
from time import time
from numpy import sqrt

start = time()
def phi_inv(q):
    return sqrt((q**2+sqrt(1+q**2))/2.0)**(-1)
def phi_inverse(q):
    return sqrt(q)*sqrt((1+sqrt(1+(4/q**2)))/2.0)
#Global parameter
a_0 = 10e-2
m = 10e11 #[MO]
R = np.arange(0.3,100,0.001) #[kpc]
G = 13.34*10e-11 #[kpc^3 MO^-1 gy^-2]
#Newton Dynamic
omega = sqrt((G*m)/(4*R**3)) #velocities
vt = omega*R
plt.xlabel("R [kpc]")
plt.ylabel("v(t)")
plt.plot(R,vt, 'r-', label="Newtonian Dynamics")
#MoND
omega = sqrt((a_0/R)*phi_inverse((m*G)/(4.0*a_0*R**2))) #velocities
vt = omega*R
plt.plot(R,vt, 'b-', label="Modified Newtonian Dynamics")
#plt.title("Total Energy of Symplectic Newton Dynamics")
#filename='./figures/plot.png'
#plt.savefig(filename)
print("Time for running ", time()-start, "seconds")
plt.legend()
plt.show()