# -*- coding: utf-8 -*-
"""
Created on Fri Feb  1 13:00:18 2019
@author: sweetc1

In Class Fri 02-10-19
"""

import numpy as np
L = 10
x = 5

def H_ij(i,j):
    a = ((np.pi**2) * (j**2)) / (2 * L**2) # kinetic energy
    b = (np.sqrt(2/L) * (np.sin((j * x * np.pi) / L))) * (np.sqrt(2/L) * (np.sin((i * x * np.pi) / L)))
    ham_int = a + b
    
    if i==j:    # include kinetic
        ham_int = a + b
    else:       # exclude kinetic
        ham_int = b
    return ham_int # with or without kinetic, depending on i & j

H_mat = np.zeros((3,3)) # 3x3 array of values

for i in range(1,4):    # fill array with values between 1 & 3
    for j in range(1,4):
        H_mat[i-1][j-1] = H_ij(i,j)
        
#print("Hamiltonian matrix:",H_mat)

# Define vectors of coefficients (c):

c = np.zeros(3) # create array 'c' to hold coefficients
c[2] = 1
norm = np.dot(np.transpose(c),c)    # norm = 1 unless coefficient product is not 1
Hc = np.dot(H_mat,c)    # product of H_mat on c
E = np.dot(np.transpose(c),Hc) / norm  # find expectation values
# transpose rotates array for multiplication purpose
#print("Energy expectation value:",E)
c_t = np.transpose(c)

print("E:",E)
#print("c_t:",c_t)
E_opt, c_opt = np.linalg.eig(H_mat)
#print("Eigenvalues:",E_opt[0])
#print("Eigenvectors:",c_opt[0])
