import matplotlib.pyplot as plt
import numpy as np
from time import time

start = time()
#Global parameter
a_0 = 10e-8
m = 10e11 #[MO]
R = np.arange(0.5,2.9,0.001) #[kpc]
G = 13.34*10e-11 #[kpc^3 MO^-1 gy^-2]
#Newton Dynamic
omega = np.sqrt((G*m)/(4*R**3)) #velocities
vt = omega*R
plt.xlabel("R [kpc]")
plt.ylabel("omega")
plt.plot(R,omega, 'r-')
#MoND
omega = np.sqrt((a_0*G*m)/(4))/R #velocities
vt = omega*R
plt.plot(R,omega, 'b-')
#plt.title("Total Energy of Symplectic Newton Dynamics")
#filename='./figures/plot.png'
#plt.savefig(filename)
print("Time for running ", time()-start, "seconds")
plt.show()