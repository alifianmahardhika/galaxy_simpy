import numpy as np
import matplotlib.pyplot as plt

#initil-parameter-symple-model
R = 2
N = 150
m = 10
m_0 = 10*m
dt = 0.002
G = 13.37*10**(-11)
eps =0.3
alpha_0 = 1.2e-8
omega = np.random.normal(0,2*np.pi)

#inital_array
x = np.zeros((N,2))
x_center = np.ones((N,2))
v = np.zeros((N,2))
f = np.zeros((N,2))

#initialization_particle-system
for i in range(N):
    t = i*dt
    x[i] = [np.random.uniform(-R,R)*np.cos(2*np.pi*(i/N)),np.random.uniform(-R,R)*np.sin(2*np.pi*(i/N))]
    v[i] = [-omega*x[i,1],omega*x[i,0]]

#force
def force(f,x,x_center,N):
    ff = np.zeros((N,2))
    for i in range(N):
        x_2nd = m_0*G*(x_center[i])
        for j in range(N):
            if i!=j:
                r_ij = 2*R*np.sin(np.pi*np.abs(x[i]-x[j])/N)
                ff = m*G*(x[i]-x[j])/(r_ij+eps)**3 - x_2nd
            f+=ff
    return f

#euler-ND
for k in range(100):
    no = np.str(k)
    x+=v*dt
    v+=force(f,x,x_center,N)*dt
    plt.xlim(-6,6)
    plt.ylim(-6,6)
    plt.scatter(x[:,0],x[:,1],s=10,c='b')
    plt.savefig('./output_euler/plot2d-euler'+no+'.png')
    plt.close()

#symplectic-ND
for k in range(100):
    no = np.str(k)
    v+=force(f,x,x_center,N)*dt
    x+=v*dt
    plt.xlim(-6,6)
    plt.ylim(-6,6)
    plt.scatter(x[:,0],x[:,1],s=10,c='r')
    plt.savefig('./output_symplectic/plot2d-symplectic'+no+'.png')
    plt.close()
    
def beta(q):  #for-MOND-model
    return np.sqrt(1+np.sqrt(1+(4/q**2))/2)
    
# symplectic_mond
a = np.zeros((N,2))
for k in range(100):
    no= np.str(k)
    qq = (np.abs(force(f,x,x_center,N))/alpha_0)
    a = force(f,x,x_center,N)*beta(qq)
    v += a*dt
    x += v*dt
    plt.xlim(-6,6)
    plt.ylim(-6,6)
    plt.scatter(x[:,0],x[:,1],s=10,c='m')
    plt.savefig('./output_mond_symplectic/plot2dmond-symplectic'+no+'.png')
    plt.close()
    
# euler_mond
a = np.zeros((N,2))
for k in range(100):
    no= np.str(k)
    qq = (np.abs(force(f,x,x_center,N))/alpha_0)
    a= force(f,x,x_center,N)*beta(qq)
    x += v*dt
    v += a*dt
    plt.xlim(-6,6)
    plt.ylim(-6,6)
    plt.scatter(x[:,0],x[:,1],s=10,c='k')
    plt.savefig('./output_mond_euler/plot2dmond-euler'+no+'.png')
    plt.close()