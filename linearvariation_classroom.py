# -*- coding: utf-8 -*-
"""
Created on Sat Jan 26 12:14:30 2019
@author: csweet

Linear Variational Method 
Full Problem Set Due Feb 15th
"""
import numpy as np
### General expression for integrals H_ij:
j = 2.
i = 1.
L = 10.
x = 5.
### Does not include d_ij conditions
psi_i = ((np.pi**2)*(j**2))/(2*(L**2))
psi_j = (1/5)*np.sin((j*x*np.pi)/L)*np.sin((i*x*np.pi)/L)
### Above are the 2 terms for the full equation below
H_ij = psi_i + psi_j

print(H_ij)
### How do I involve the loop to dictate d_ij???
"""
### use this loop to determine d_ij:
def H_ij(i,j):
    
    if i==j:
    ### if i=j, d_ij=1
    ### include first term of equation    
    else:
    ### if i=/=j, d_ij=0
    ### exclude first term    
    return ham_int
"""
