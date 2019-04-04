#!/usr/bin/env python
# run in a DA dir (has files: xas.bin, xbs.bin, obs.csv),
# QuickPlots.py <truth.bin>
import csv
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

def help(): #{{{
    print('''
QuickPlotATimeStep.py <truth.bin>. <t, the time step index to plot, default 0>, <nx, default 1000> (Note, run in slam or 4dvar running dir)
''')
    exit(0)
#}}}


# return dict: obs_t -> [ (obs_x, obs_value, obs_error) ] 
def read_obs_dict(obs_file):  #{{{
    dct={}
    with open(obs_file, "r") as f:
       reader = csv.DictReader(f)
       for row in reader:
            t=int(row["t"])
            if(t not in dct):
                dct[t] = []
            x=int(row["x"])
            val = float(row["obs"])
            err = float(row["err_stdv"])
            dct[t].append((x,val,err))
    return dct
#}}}

# all keys of obs_dict, minus t
def offset_obs_dict(obs_dict, t):
    new_dct = {}
    for key in obs_dict.keys():
        new_dct[key - t] = obs_dict[key]
    return new_dct


#full field of Xb, Xa, Xt, Obs for each time
#problem is full field has much larger magnitude than err
def plot_atime_xbxaxt_obs(xts, xbs, xas, obs_dict, t_step=1): #{{{
    nt,nx = xts.shape
    for t in range(0, nt, t_step):
        if(t not in obs_dict):
            continue
        print("--t=%d"%t)
        fig=plt.figure(figsize=(10,6))
        xcoord = np.arange(0,nx)
        plt.plot(xcoord, xts[t,:], 'k--', label="xt")
        plt.plot(xcoord, xbs[t,:], 'b', label='xb')
        plt.plot(xcoord, xas[t,:], 'r', label='xa')
        obs_list = obs_dict[t]  #[(obs_x, obs_val, obs_err) ... ]
        obs_xs, obs_vals, obs_errs = zip(*obs_list)
        plt.plot(obs_xs, obs_vals, 'g*', label='obs')
        plt.legend()
        plt.xlabel('x coord')
        plt.savefig("images/xbxaxt_obs_t%0.2d"%t)
        plt.close()
#}}}

#similar to plot_atime_xbxaxt_obs, but plot the error (xb/xa/obs - xt), to emphasize the error
def plot_atime_xbxaobs_err(xts, xbs, xas, obs_dict, t_step=1):
    nt,nx = xts.shape
    for t in range(0, nt, t_step):
        if(t not in obs_dict):
            continue
        print("--t=%d"%t)
        fig=plt.figure(figsize=(10,6))
        xcoord = np.arange(0,nx)
        zeros = np.zeros((nx,))
        plt.plot(xcoord, zeros, 'k--')
        plt.plot(xcoord, xbs[t,:] - xts[t,:], 'b', label='xb-xt')
        plt.plot(xcoord, xas[t,:] - xts[t,:], 'r', label='xa-xt')
        obs_list = obs_dict[t]  #[(obs_x, obs_val, obs_err) ... ]
        obs_xs, obs_vals, obs_errs = zip(*obs_list)
        obs_mt = np.array(obs_vals) - xts[t,np.array(obs_xs)]
        plt.plot(obs_xs, obs_mt, 'g*', label='obs-H(xt)')
        plt.legend()
        plt.xlabel('x coord')
        plt.savefig("images/xbxaobs_err_t%0.2d"%t)
        plt.close()


def check_file_exist(files): #{{{
   for filex in files:
       if(not os.path.isfile(filex)):
           print("Error, %s doesn't exist!" %filex)
           exit(1)
#}}}
    


if __name__ == "__main__":
    n_arg = len(sys.argv[1:])
    nx=0
    it=0
    if(n_arg == 1):
        xts_file = sys.argv[1]
        nx=1000
        it=0
    elif(n_arg == 3):
        xts_file = sys.argv[1]
        it=int(sys.argv[2])
        nx = int(sys.argv[3])
    elif(n_arg == 2):
        xts_file = sys.argv[1]
        it=int(sys.argv[2])
        nx = 1000
    else:
        help()

    check_file_exist(["obs.csv", "xbs.bin", "xas.bin"])
    
    print("start reading files")
    xbs=np.fromfile("xbs.bin", dtype="double", count=-1)
    n = xbs.shape[0]
    nt = int(n/nx)
    xbs=xbs.reshape((nt,nx))
    xb0 = xbs[it,:].reshape(1,nx)
    
    xas=np.fromfile("xas.bin", dtype="double", count=-1)
    nt = int(xas.shape[0] / nx)
    xas=xas.reshape((nt,nx))
    xa0 = xas[it,:].reshape(1,nx)

    xts=np.fromfile(xts_file, dtype="double", count=-1)
    nt = int(xts.shape[0] / nx)
    xts=xts.reshape((nt,nx))
    xt0 = xts[it,:].reshape(1,nx)

    obs_dict = offset_obs_dict(read_obs_dict("obs.csv"), it)

    if(not os.path.isdir("images")):
        os.mkdir("images")


    print("start plotting Xb/Xa/Xt/Obs")
    plot_atime_xbxaxt_obs(xt0, xb0, xa0, obs_dict, 36) #36=6*6, every 6 hour one plot
    plot_atime_xbxaobs_err(xt0, xb0, xa0, obs_dict, 36)

    

