#!/usr/bin/env python
import sys
import numpy as np
from scipy import integrate

def burgers_spatialcentral(t, u, dx, R):
    nx = len(u)
    fu = np.zeros((nx,), np.double)
    for i in range(nx):
        ip1 = (i+1)%nx
        im1 = (i-1+nx)%nx
        fu[i] = -u[i]*(u[ip1] - u[im1])/(2*dx) + (u[ip1]+u[im1]-2*u[i]) / (R*dx*dx)
    return fu

def help():
   print('''
python burgers.py [method: lsoda, dopri5] [init file] [nt](will output 0..nt, total nt+1 steps)
''')
   exit(0)


if __name__ == "__main__":

    if (len(sys.argv[1:]) < 3):
        help()

    method=sys.argv[1]
    init_file=sys.argv[2]
    nt=int(sys.argv[3])

    x0=np.fromfile(init_file, dtype="double", count=-1)
    
    dt=600.
    dx=40000.
    R=1e-6
    t=np.arange(0, dt*nt+1, dt)
    #lsoda
    #r = integrate.ode(burgers_spatialcentral).set_integrator('lsoda', first_step=dt)
    #dopri5
    r = integrate.ode(burgers_spatialcentral).set_integrator('dopri5', first_step=dt)
    r.set_initial_value(x0).set_f_params(dx,R)
    
    xs=[x0]
    cnt=0
    while r.successful() and r.t < nt*dt:
        cnt+=1
        xs.append(r.integrate(r.t+dt))
    
    xs=np.array(xs)
    xs.tofile("output.bin")


