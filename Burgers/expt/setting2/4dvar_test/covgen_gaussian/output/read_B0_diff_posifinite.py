#!/usr/bin/env python

import numpy as np

BCov=np.fromfile("Cov_B0.bin", dtype="double", count=-1)
lenB = len(BCov)
nx = int(np.sqrt(lenB))
Binflator=0.04
BCov = np.reshape(BCov,(nx,nx))*Binflator

def is_matrix_positve_finite(A):
    try:
        Z=np.linalg.cholesky(A)
        return True
    except:
        return False
print(BCov[0:5,0:5])

# increase lambda from 1, to make Amatrix = lambda * diag(BCov) - BCov, positive finite
lambdas=np.arange(120,130,1)
for lamb in lambdas:
    Amtx = lamb * np.diag(np.diag(BCov)) - BCov
    is_posfin = is_matrix_positve_finite(Amtx)
    print("lambda = %d" %lamb + str(is_posfin))
    if (is_posfin):
        print(Amtx[0:5,0:5])
        break



