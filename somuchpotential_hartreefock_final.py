"""
So Much Potential
02.15.19
"""

import numpy as np
from matplotlib import pyplot as plt
from scipy.interpolate import InterpolatedUnivariateSpline

### create an array of 20 bond lengths spanning 0.5 - 3.5 angstroms
### but store in atomic units of length... note 1 angstrom = 1.88973 a.u. of length
r_array = np.linspace(0.5,3.5,20)*1.88973
print(r_array)

### array below is calculated in Jupyter notebook
### these values are fed into file on Janet
### resulting values are added to array underneath
"""
[0.944865   1.24324342 1.54162184 1.84000026 2.13837868 2.43675711
 2.73513553 3.03351395 3.33189237 3.63027079 3.92864921 4.22702763
 4.52540605 4.82378447 5.12216289 5.42054132 5.71891974 6.01729816
 6.31567658 6.614055  ]
"""
### fill in this array with your actual energy values... there should be 20 values in total!!!
### these bond length values are taken from the output files generated on Janet 
E_array = [-107.141,-110.874, -112.204, -112.622, -112.699, -112.653, -112.570, 
          -112.484, -112.410, -112.361, -112.333, -112.315, -112.304, -112.069,
          -112.291, -112.287, -112.001, -111.948, -111.927, -112.280]
# E_array is array of bond lengths

### "wiggle" in the graph is consequence of imperfect approx. by H-F
# TO see the plot of the PES, uncomment the following lines
plt.plot(r_array, E_array, 'red')
plt.show()

### use cubic spline interpolation
order = 3
### form the interpolator object for the data
sE = InterpolatedUnivariateSpline(r_array, E_array, k=order)
### form a much finer grid
r_fine = np.linspace(1.06,5.0,200)
### compute the interpolated/extrapolated values for E on this grid
E_fine = sE(r_fine)
### plot the interpolated data
plt.plot(r_fine, E_fine, 'blue')
plt.show()

### take the derivative of potential
fE = sE.derivative()
### force is the negative of the derivative
F_fine = -1*fE(r_fine)

### plot the forces
plt.plot(r_fine, np.abs(F_fine), 'black')
plt.xlim(1,5)
plt.show()

### Find index of the PES where it has its global minimum
r_eq_idx = np.argmin(E_fine)
### find the value of the separation corresponding to that index
r_eq = r_fine[r_eq_idx]
### print equilibrium bond-length in atomic units and in angstroms
print("Equilibrium bond length is ",r_eq," atomic units")
print("Equilibrium bond length is ",r_eq*0.529," Angstroms")

### get second derivative of potential energy curve... recall that we fit a spline to
### to the first derivative already and called that spline function fE.
cE = fE.derivative()

### evaluate the second derivative at r_eq to get k
k = cE(r_eq)
"""
Q1.0.1:
    What is the reduced mass of the CO milecule in atomic units?
"""
# Proton = 1836 amu
C_mass = 12 * 1836   # mass of carbon atom in amu
O_mass = 16 * 1836   # mass of oxygen atom in amu
CO_mu = (C_mass * O_mass) / (C_mass + O_mass) # calculate reduced mass
print("Calculated reduced mass of CO molecule:",CO_mu)
### Reduced mass of CO molecules is 6.857 amu

### define reduced mass of CO as m_C * m_O /(m_C + m_O) where mass is in atomic units (electron mass = 1)
mu = 13625. # answer in note above is different...

### define "ground-state" velocity:
v = np.sqrt( np.sqrt(k/mu)/(2*mu))

### get random position and velocity for CO within a reasonable range
r_init = np.random.uniform(0.75*r_eq,2*r_eq)
v_init = np.random.uniform(-2*v,2*v)

### print initial position and velocity
print("Initial separation is ",r_init, "atomic units")
print("Initial velocity is   ",v_init, "atomic units")
### establish time-step for integration to be 0.2 atomic units... this is about 0.01 femtoseconds
dt = 0.02

### get force on particle 
F_init = -1*fE(r_init)

def Velocity_Verlet(r_curr, v_curr, mu, f_interp, dt):
    ### get acceleration at current time
    a_curr = -1*f_interp(r_curr)/mu
    
    ### use current acceleration and velocity to update position
    r_fut = r_curr + v_curr * dt + 0.5 * a_curr * dt**2
    
    ### use r_fut to get future acceleration a_fut
    ''' STUDENT WRITTEN CODE GOES HERE!'''
    ### look at jupyter notebook for equation
    a_fut = -1 * f_interp (r_fut) / mu
    ### use current and future acceleration to get future velocity v_fut
    ''' STUDENT WRITTEN CODE GOES HERE!'''
    ### look at jupyter notebook for equation
    v_fut = v_curr + 0.5 * (a_curr + a_fut) * dt
    
    result = [r_fut, v_fut]
    
    return result

"""
Q1.0.2:
    Use spline fit of the PES of the OC molecule to estimate the vibrational
    frequency of CO. Express your number in atomic units and also convert to a
    common spectroscopic unit system of your choosing (wavenumbers, nm, um,
    THz, are al acceptable choices)
"""    
r = 2 * np.pi * np.sqrt(CO_mu / k)
print("Vibrational frequency (atomic units):",r)
print("Vibrational frequency (nm):",r * 0.053)


### how many updates do you want to perform?
N_updates = 200000

### Now use r_init and v_init and run velocity verlet update N_updates times, plot results
### these arrays will store the time, the position vs time, and the velocity vs time
r_vs_t = np.zeros(N_updates)
v_vs_t = np.zeros(N_updates)
t_array = np.zeros(N_updates)

### first entry is the intial position and velocity
r_vs_t[0] = r_init
v_vs_t[0] = v_init

### first Velocity Verlet update
result_array = Velocity_Verlet(r_init, v_init, mu, fE, dt)

### do the update N_update-1 more times
for i in range(1,N_updates):
    tmp = Velocity_Verlet(result_array[0], result_array[1], mu, fE, dt)
    result_array = tmp
    t_array[i] = dt*i
    r_vs_t[i] = result_array[0]
    v_vs_t[i] = result_array[1]

### Plot the trajectory of bondlength vs time:
plt.plot(t_array, r_vs_t, 'red')
plt.show()

### plot the phase space trajectory of position vs momentum
plt.plot(mu*v_vs_t, r_vs_t, 'blue')
plt.show()

"""
Q1.0.3:
    What will be the acceleration of the bond stretch when C is separted from O
    by 3 atomic units? You can express your acceleration in atomic units, also.
"""
a = -1 * fE(3) / CO_mu 
print("Bond acceleration at 3 atomic units:",a)