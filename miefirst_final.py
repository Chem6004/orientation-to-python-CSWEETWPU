"""
Mie First
12-Apr-19
"""
"""
Silver: dE = E_2 - E_1 ~ 3.5 eV
mu_12 = mu_21 ---> energy of excitement is same as energy of de-excitement

H_t = H_o - E_t * mu ---> mu is a potential
***When light is on, potential is time-dependent***

mu_t = Tr(D_t * mu)
D = psi * psi_star ---> density function
psi ---> vector of coefficients (c_1, c_2, etc.)
D = [c_1 * c_1star, c_1 * c_2star
     c_2 * c_1, c_2 * c_2star]
***Matrix corresponds to coefficients of wavefunction***

If psi = 1 * psi_1_x + 0 * psi_2_x:
Matrix = [0, 0
          0, 1]

If psi = sqrt(1/2) * psi_1 + sqrt(1/2) * psi_2:
Matrix = [1/2, 1/2,
          1/2, 1/2]
***Diagonals must equal to 1***

TDSE for density matrix = Liouville equation:
    d/dt D = -1/h_bar (HD-DH)
***Order of D & H do matter ---> commutater***
    
If H_ij = <psi_i|H_hat|psi_j>:
Matrix = [E_1,0,
          0, E_2]

E_1 = 0 (ground state energy), E_2 = excited state energy
***The relative gap between E_1 & E_2 is what is important***        

***If d/.dt D = 0, matrice commute = no time derivative***
***If d/dt D =/= 0, do not commute = there is a time derivative***

Apply Euler method:
D(t+dt) = D_t + d/dt D_t * dt
"""

from matplotlib import pyplot as plt
import numpy as np
from numpy import linalg as LA
import math


def Euler(H0, mu, Vint, gamma, D, h, t, tau):
### H0 is matrix passesd to the fx, EField is matrix that comes out
### Defined Hamiltonian at current time (EField changes with time)    
    H = H0 - EField(t, tau) * mu
    ### Defined time derivative of density matrix at current time by evaluating
    ### Liouville-Linblad equation
    Ddot = Liouville(H, D) + Lindblad(D, gamma)
    ### New density matrix defined as density matrix at current time + time
    ### derivative multiplied by timestep (h) based on Euler update
    Dnew = D + h*Ddot
    
    return Dnew


def Lindblad(D, gamma):
    dim = len(D)
    LD = np.zeros_like(D)
    ### need |g><g|
    bra_1 = CreateBas(dim, 0)
    gm = Form_Rho(bra_1)
    
    for k in range(1,dim):
        bra_k = CreateBas(dim, k)
        km = Form_Rho(bra_k)
        
        ### first term 2*gam*<k|D|k>|g><g|
        t1 = 2*gamma*D[k][k]*gm
        ### second term is |k><k|*D
        t2 = np.dot(km,D)
        ### third term is  D*|k><k|
        t3 = np.dot(D, km)
        LD = LD + t1 - gamma*t2 - gamma*t3
        
    return LD

### Take commutator of H and D to give Ddot
def Liouville(H, D):
    ci = 0.+1j
    return -ci*(np.dot(H,D) - np.dot(D, H))

def EField(t, tau):
    Ef = 0.
    if t<tau:
        Ef = 0.001*np.sin(t*np.pi/tau)*np.sin(t*np.pi/tau)*np.sin(0.1192*t)
    return Ef

def Form_Rho(Psi):

    D = np.outer(Psi,np.conj(Psi))
    return D

### Creates basis vector for state k
### k=0 -> ground state, k=1 -> first excited-state, etc
def CreateBas(dim, k):
    bas = np.zeros(dim)
    bas[k] = 1
    return bas

### Runge-Katta Method

def RK4(H0, mu, Vint, gamma, D, h, t, tau):
    k1 = np.zeros_like(D)
    k2 = np.zeros_like(D)
    k3 = np.zeros_like(D)
    k4 = np.zeros_like(D)
    D1 = np.zeros_like(D)
    D2 = np.zeros_like(D)
    D3 = np.zeros_like(D)
    D4 = np.zeros_like(D)
    Df = np.zeros_like(D)
    
    ### Get k1
    H1 = H0 - EField(t, tau)*mu + Vint
    D1 = D    
    k1 = h*Liouville(H1,D1) + h*Lindblad(D1, gamma)
    
    ## Update H and D and get k2
    H2 = H0 - EField(t+h/2, tau)*mu + Vint
    D2 = D+k1/2.
    k2 = h*Liouville(H2, D2) + h*Lindblad(D2, gamma)
    
    ### UPdate H and D and get k3
    H3 = H2
    D3 = D+k2/2
    k3 = h*Liouville(H3, D3) + h*Lindblad(D3, gamma) 
    
    ### Update H and D and get K4
    H4 = H0 - EField(t+h, tau)*mu + Vint
    D4 = D+k3
    k4 = h*Liouville(H4, D4) + h*Lindblad(D4, gamma)
    
    Df = D + (1/6.)*(k1 + 2.*k2 + 2*k3 + k4)
    return Df

### Third piece

### Set up some parameters for the quantum dynamics simulation
dt = 0.1
tau = 100 #150.
gamma = 0.0017
eps0 = 8.854e-12
mu_au_to_si = 8.47835326e-30
E_au_to_si = 5.14220652e11
mu_z = 58.

### Create some arrays
MUZ= np.zeros((2,2),dtype=complex)
Vint = np.zeros((2,2),dtype=complex)
### Density matrix for RK4 updates
D_RK4 = np.zeros((2,2),dtype=complex)
### Density matrix for Euler updates
D_EU  = np.zeros((2,2),dtype=complex)
H0 = np.zeros((2,2))

### initialize values of the arrays for Hamiltonian and Density matrices
H0[0][0] = 0.1275
D_RK4[0][0] = 1.+0j
D_EU[0][0] = 1.+0j
MUZ[0][1] = mu_z
MUZ[1][0] = mu_z

### create arrays for time-dependent quantities
Nsteps = 40000
ez = np.zeros(Nsteps)
### array for mu(t) for RK4 updates
mu_of_t_rk4 = np.zeros(Nsteps,dtype=complex)
### array for mu(t) for Euler updates
mu_of_t_eu = np.zeros(Nsteps,dtype=complex)
time = np.zeros(Nsteps)
energy = np.zeros(Nsteps)


### Run the dynamics
for i in range(0,Nsteps):
    ### reciprocal axis
    energy[i] = np.pi*2*(i+1)/(Nsteps*dt)
    ### time access
    time[i] = i*dt
    ### time-dependent electric field
    ez[i] = EField(i*dt, tau)*E_au_to_si
    ### update to the Density matrix using RK4
    D_RK4 = RK4(H0, MUZ, Vint, gamma, D_RK4, dt, dt*i, tau)
    D_EU = Euler(H0, MUZ, Vint, gamma, D_EU, dt, dt*i, tau)
    ### Update to mu(t) using RK4
    DMU_RK4 = np.matmul(D_RK4, MUZ)
    DMU_EU = np.matmul(D_EU, MUZ)
    mu_of_t_rk4[i] = (DMU_RK4[0][0] + DMU_RK4[1][1])*mu_au_to_si
    mu_of_t_eu[i] = (DMU_EU[0][0] + DMU_EU[1][1])*mu_au_to_si
    ### add update using Euler step!!!

print(energy[0])

plt.plot(time, mu_of_t_rk4, 'red', time, mu_of_t_eu, 'b--', time, ez, 'green')
plt.show()

"""
mu_freq_rk4 = np.fft.fft(mu_of_t_rk4)/(Nsteps)
ez_freq = np.fft.fft(ez)/(Nsteps)
alpha_rk4 = mu_freq_rk4/ez_freq
lam = 1e-9*1240/(energy*27.211) ### in nm

sigma_rk4 = 2*np.pi/(lam*eps0) * np.imag(alpha_rk4)
plt.plot(energy*27.211, sigma_rk4, 'red')
#plt.plot(1240/(sphere.lambda_array*1e9), sigma_abs, 'b--')
plt.xlim(1.5,4.0)
plt.ylim(0,1e-16)
#plt.plot(time, ez, 'red', time, mu_of_t, 'blue')
"""
plt.show()