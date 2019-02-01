# -*- coding: utf-8 -*-
"""
Created on Sat Jan 26 12:39:08 2019
@author: csweet

Orientation to Python Problems
Due Fri Feb 08th
"""
### Import libraries needed later on up here:
import numpy as np
import time
from matplotlib import pyplot as plt

"""
Question #1:
Hows does total KINETIC energy of a collection of N particles
grow with N assuming each particle has the same average kinetic
energy? Compute the total kinetic energy for 5 different values
of N & polt the results using 'pyplot'.
"""

### Initialize variables
N = 5               # N, number of particles
m = np.zeros(N)     # Array of mass (zero filled) of size N
v = np.zeros(N)     # Array of velocities (zero filled) of size N

### Assign values to mass and velocity arrays
for i in range(0,N):
    m[i] = 1        # m_i = 1
    v[i] = 2.5      # v_i = 2.5
    
### Un-comment to print the arrays defined above:
#print("Array of masses: ",m)
#print("Array of velocities: ",v)

### COMPUTE kinetic energy of particles
T = 1/2 * m * v**2      # T is an array bc m & v ar arrays
#print("Array of kinetic energy:",T)               # Print T to prove it is an array

### SUM kinetic energy of particles
"""
T_loop = 0              # Can use nested loops (convoluted)

for i in range(0,N):
    T_loop = T_loop + T[i]
"""
    
#T_sum = np.sum(T)       # Use of np.sum function is easier

#print("Loop sum:",T_loop)
#print("Numpy sum:",T_sum)

### PLOT kinetic energy of particles
### x-values are N, y-values are resulting T
S = np.zeros(N)
for i in range(0,N):
    if i==0:
        S[i] = T[i]
    else:
        S[i] = T[i] + S[i-1]
        
#print("Discrete kinetic energies:",S)
#print("\nKinetic Energy (red):")
#plt.plot(np.linspace(0,1,N),S,'red')
### Hooray, its a linear array!!!

### Define charge (q), position (x), & separation (r)
q = np.ones(N)                  # Array filled with all 1's
x = np.linspace(0,(N-1)*0.2,N)  # Start, increment, end
r = np.zeros((N,N))             # 2-Dimensional array N x N

### Compute separations using nested for loops
for i in range(0,N):            # First particle interacting
    for j in range(0,N):        # Second particle interacting
        r[i][j] = np.sqrt((x[i]-x[j])**2)
        
#print("Array of charges (q):",q)
#print("Array of positions (x):",x)
#print("Array of separations (r):",r)

"""
Question #2:
How does the total POTENTIAL energy of a collection of N
equally-spaced charged particels grow with N? COmput the total
potential energy for 5 different values on N & plot the results.
"""
### Use array of separations to compute potential energy
### Create a function called potential to pass r & q to
def potential(sep_array,chg_array,N):
    pot = 0
    # Use nested loops
    for i in range(0,N):
        for j in range(0,N):
            if (i!=j):
                pot = pot + chg_array[i] * chg_array[j] / sep_array[i][j]
    return pot

S_v = np.zeros(N)
for i in range(0,N):
    V = potential(r,q,i)
    if i==0:
        S_v[i] = V
    else:
        S_v[i] = V + S_v[i-1]
    

#print("Discrete potentials:",S_v)
#print("Potential Energy (blue):")
#plt.plot(np.linspace(0,1,N),S_v,'blue')

"""
Question 3:
Use the 'time' library to determine the time required to
compute the KINETIC & POTENTIAL energy for the 5 different
values of N; plot the time required vs. N & discuss if the
kinetic energy seems to scale linearly & the potential energy
seems to scale quadratically with N.
"""
start = time.time()
x = np.linspace(1,5,5)
y = S
z = S_v

plt.plot(x,y,'red')
print("\nKinetic energy over time (red)")
plt.plot(x,z,'blue')
print("Potential energy over time (blue)")

end = time.time()
l = end - start
print("Total run time:",l,"sec")

"""
Example pyplot using time: apply this to exercises above
   1 - Create array of x-values (start pos) & y-values (end pos)
   2 - Measure the time it takes to run the entire program
"""
"""
start = time.time()             # Get time at beginning
x = np.linspace(-5,5,100)       # 100 x-values btw -5 & 5
y = x**2                        # Array of y-values

#plt.plot(x,y,'red')             # Graphs x vs y in red
plt.show()

end = time.time()
length = end - start

#print("Total run time (sec) is:",length)
"""