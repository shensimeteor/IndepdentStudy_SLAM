#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
file="../burgers_run.bin"
nx=1000
N=101
Xs=np.fromfile(file, dtype="double", count=-1).reshape(N,nx)

#print(Xs)
plot_dt=10

xcoord=np.arange(0,nx)
for i in range(0,N,plot_dt):
    print(i)
    plt.figure()
    plt.plot(xcoord, Xs[i][:], '-o')



