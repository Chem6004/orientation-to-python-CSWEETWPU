# -*- coding: utf-8 -*-
"""
Created on Thu Jan 31 22:52:14 2019
@author: csweet

Linear Variational Method
Problems Due Fri Feb 12, 2019
"""
import numpy as np
"""
Pt1 Q1:
Work out general expression for integral H_ij
"""
N = 3                       # Number of particles
L = 10                      # length of 'box' for particle in a box
x = 5                       # delta(x-5) = sub x for 5
i = np.zeros(N)             # array for first particle i
j = np.zeros(N)             # array for second particle j

for n in range(0,N):        # should increase i & j by 1 for N times
    i[n] = 1 
    j[n] = 1

print("array of i:",i)
print("array of j:",j)

i_d = np.zeros(N)           # assign discrete values to i
for n in range(0,2):
    if n==0:
        i_d[n] = i[n]
    else:
        i_d[n] = i[n] + i_d[n-1]

j_d = np.zeros(N)           # assign discrete values to j
for n in range(0,2):
    if n == 0:
        j_d[n] = j[n]
    else:
        j_d[n] = j[n] + j_d[n-1]
        
print("array of discrete i:",i_d)
print("array of discrete j:",j_d)

fst = ((np.pi**2)*(j**2))/(2*(L**2))
sec = ((1/5)*np.sin((j*x*np.pi)/L))*(np.sin((i*x*np.pi)/L))

d_ij = np.zeros(N)          # decide if d_ij is 1 or 0
for n in range(0,1):
    if n == n:
        d_ij[n] = 1
    else:
        d_ij[n] = 0
        
H_ij = (fst*d_ij) + sec     # compute integral value

print("array of integrals:",H_ij)
print("array of first terms:",fst*d_ij)
print("array of second terms:",sec)
