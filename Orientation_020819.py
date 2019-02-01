# -*- coding: utf-8 -*-
"""
Created on Mon Jan 28 21:25:11 2019

@author: csweet
"""

### Orientation to Python
### Due Fri Feb 08, 2019

import numpy as np
import time
from matplotlib import pyplot as plt

# Number of particles:
N = 5

### Question #1:
"""
Hows does total KINETIC energy of a collection of N particles
grow with N assuming each particle has the same average kinetic
energy? Compute the total kinetic energy for 5 different values
of N & polt the results using 'pyplot'.
"""

# Array of masses:
m = np.zeros(N)
# Array of velocities:
v = np.zeros(N)

# Assign values to mass & velocity arrays:
for i in range(0,N):
    m[i] = 1
    v[i] = 2.5
    
# Compute kinetic energies of particles:
K = (1/2)*m*(v**2)

# Sum kinetic energies:
K_sum = np.sum(K)

# Get discrete kinetic energy values to plot:
K_dis = np.zeros(N)

for i in range(0,N):
    if i==0:
        K_dis[i] = K[i]
    else:
        K_dis[i] = K[i] + K_dis[i-1]
  
      
print("\n")
print("Masses:",m)
print("Velcoities:",v)
print("Kinetic Energies:",K)
print("Total Kinetic Energy:",K_sum)
print("Discrete Kinetic Energies:",K_dis)

### Plot kinetic energy vs. number of particles:
plt.plot(np.linspace(0,1,N),K_dis,'red')
print("\nRed = Kinetic Energy vs. Particles")

### Question #2:
"""
How does the total POTENTIAL energy of a collection of N
equally-spaced charged particels grow with N? COmput the total
potential energy for 5 different values on N & plot the results.
"""

