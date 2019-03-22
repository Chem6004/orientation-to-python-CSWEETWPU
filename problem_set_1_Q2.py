"""
Chelsea Sweet
Problem Set #1
Question #2
Dr. Foley
"""

import numpy as np
x = 1

"""
Q-1-1:
Work out a general expression for integras H_ij
"""
def H_ij(i,j):
    a = (np.pi**2) * (j**2) # kinetic energy
    b = x * (x - 1)
    ham_int = a + b
    
    """
    Q-1-2:
    Write a py fxn that takes the indices i & j and returns the value of H_ij
    """    
    if i==j:    # include kinetic
        ham_int = a + b
    else:       # exclude kinetic
        ham_int = b
    return ham_int # with or without kinetic, depending on i & j

H_mat = np.zeros((3,3)) # 3x3 array of values

for i in range(1,4):    # fill array with values between 1 & 3
    for j in range(1,4):
        H_mat[i-1][j-1] = H_ij(i,j)
        
print("Hamiltonian matrix:",H_mat)

# Define vectors of coefficients (c):

c = np.zeros(3) # create array 'c' to hold coefficients
c[0] = 1 #change the number in brackets to change E
norm = np.dot(np.transpose(c),c)    # norm = 1 unless coefficient product is not 1
Hc = np.dot(H_mat,c)    # product of H_mat on c
E = np.dot(np.transpose(c),Hc) / norm  # find expectation values
# transpose rotates array for multiplication purpose
#print("Energy expectation value:",E)
c_t = np.transpose(c)

print("Hc matrix:",Hc)
print("Expectation value:",E)

