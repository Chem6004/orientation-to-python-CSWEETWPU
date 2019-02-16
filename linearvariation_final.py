"""
Created on Fri Feb  1 13:00:18 2019
@author: sweetc1
In Class Fri 02-10-19
Due Fri 02-15-19
"""

import numpy as np
L = 10
x = 5

"""
Q-1-1:
Work out a general expression for integras H_ij
"""
def H_ij(i,j):
    a = ((np.pi**2) * (j**2)) / (2 * L**2) # kinetic energy
    b = (np.sqrt(2/L) * (np.sin((j * x * np.pi) / L))) * (np.sqrt(2/L) * (np.sin((i * x * np.pi) / L)))
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

H_mat = np.zeros((6,6)) # 3x3 array of values

for i in range(1,7):    # fill array with values between 1 & 3
    for j in range(1,7):
        H_mat[i-1][j-1] = H_ij(i,j)
        
#print("Hamiltonian matrix:",H_mat)
        
"""
Q-1-3:
Show that differentiating the energy functional w.r.t. all coefficients and
setting the derivative to 0 results in the eigenvalue eq. H_c = E_c

Not really sure how to answer this question...
"""
# Define vectors of coefficients (c):

c = np.zeros(6) # create array 'c' to hold coefficients
c[0] = 1 #change the number in brackets to change E
norm = np.dot(np.transpose(c),c)    # norm = 1 unless coefficient product is not 1
Hc = np.dot(H_mat,c)    # product of H_mat on c
E = np.dot(np.transpose(c),Hc) / norm  # find expectation values
# transpose rotates array for multiplication purpose
#print("Energy expectation value:",E)
c_t = np.transpose(c)

print("Hc matrix:",Hc)
print("Expectation value:",E)


#print("E:",E)
#print("c_t:",c_t)
E_opt, c_opt = np.linalg.eig(H_mat)
#print("Eigenvalues:",E_opt[0])
#print("Eigenvectors:",c_opt[0])

"""
Q-2-1:
Is the energy calculated higher or lower that ground-state energy of the
ordinary PIB system w/o delta fxn?
"""
E_sansdelta = (np.pi**2)/(2*(L**2))

#print("Ground-state E of PIB w/o delta potential:",E_sansdelta)

### The energy calculated above is greater than the PIB gorund state
### Can't go lower than ground state, so this is the expected outcome
#print(E,"is greater than",E_sansdelta)
"""
Q-2-2:
Why do you think that mixing in fxns that correspond to excited states in the
ordinary PIB system helped improve/lower energy with the delta fxn potential?
"""
### Adding an excited state wfxn spreads the energy out so it is no longer
### concentrated at the potential spike, lowering overall energy of the system
"""
Q-2-3:
Increase the number of basis fxns to 6 (makes H a 6x6 matrix & C has 6 vector
entries) and repeat the calculation for variational estimate of ground state
energy. Does the energy improve/lower compared to with 3 basis fxns?

Hc matrix (3x3): [ 1.24674011e+00  1.22464680e-16 -1.00000000e+00]
Expectation value: 0.24934802200544673

Hc matrix (6x6): [ 1.24674011e+00  1.22464680e-16 -1.00000000e+00 -2.44929360e-16
  1.00000000e+00  3.67394040e-16]
Expectation value: 0.24934802200544673

Expectation value matches the first Hc value when c[0] = 1
"""
### The energy does not change when the matrix size is increadsed to 6x6
### The first 3 values match the 3x3 matrix Hc values