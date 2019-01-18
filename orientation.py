# -*- coding: utf-8 -*-
"""
Created on Fri Jan 18 13:39:36 2019

@author: sweetc1
"""
#print("Hello World!")

### Ignore me
"""
Ignore
Everything
"""

### These lines are class notes!
#these are commented out pieces of code

import numpy as np
### Variables for particle 1
Npart = 10
### Create empty lists for each quantity
### Zeros fxn is just number of values, not actual values
m = np.zeros(Npart)
v = np.zeros(Npart)
q = np.zeros(Npart)
x = np.zeros(Npart)
T = np.zeros(Npart)
### Print the array of zeros for m
print(m)
### Use a for loop to assign value to the variables with arrays
### Everything within the loops hould have same indentation!
### i is the variable being changed with this loop, values starting at 0
for i in range(0, Npart):
    m[i] = 1.0
    q[i] = 1.0
    x[i] = 0.5 * i
    v[i] = 1.2 * i
    ### Can use existing for loop to calculate kinetic energy
    ### Now that mass & velocity is calculates
    T[i] = 0.5 * m[i] * v[i]**2
    ### Order of calculated values within a loop matters!
print(T)
### No closing parentheses required, just chage indent   

    
### Velocity in this world is 1-D scalar
### Need at least 2 particles to find potential energy
### Need to know their charges (q) & the distance between them (r)

#T1 = 0.5 * m * v**2
### x**2 = x^2 in python language
#T2 = 1/2 * m *v * v
#T3 = 0.5 * m * v * v
#T4 = m * v**2 / 2

### Print fxn passes stored values of kinetic energy (should all be the same)
#print(T1, T2, T3, T4)
### Variables for particle 2
#m2 = 1.0
#v2 = 2.5
#q2 = 1.0
#x2 = 0.5

### Variables for pair of particles
r12 = np.sqrt((x1 - x2)**2)
V12 = (q1 * q2) / r12

print("Separation:", r12)
print("Potential Energy:",V12)
### This problem is simple enough to solve by hand...
### For many, many french fries (lol), we can make an array: list of numbers
### Array can contain mass, velocity, kinetic energy, etc. for each particle
### Take the sume of all T to get the total kinetic energy