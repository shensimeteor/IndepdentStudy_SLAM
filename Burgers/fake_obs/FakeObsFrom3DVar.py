#!/usr/bin/env python

import numpy as np
import pprint as pp
import csv
import json
import sys
import matplotlib.pyplot as plt


def help():
    print( '''
usage: FakeObsFrom3DVar.py [fakeobs_conf.json]
''')
    exit(0)

def output_csv(res, outfile):
    with open(outfile, "w") as f:
        writer = csv.writer(f)
        writer.writerow(["t","x", "obs", "err_stdv"])
        for tup in res:
            writer.writerow(tup)

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

# return 2d array xs
def read_xs_bin(xs_file, nx):
    xs=np.fromfile(xs_file, dtype="double", count=-1)
    n = xs.shape[0]
    nt = int(n/nx)
    xs=xs.reshape((nt,nx))
    return xs

# return a list of  [(it, ix, value, noise)]
def makeFakeObs(dxa0, xb0, t, xb0_err_stdv, obs_err_stdv):
    nx=xb0.shape[0]
    obs_list =  []
    for iobs in range(nx):
        obs_val = xb0[iobs] + (xb0_err_stdv**2 + obs_err_stdv**2) / (xb0_err_stdv**2) * (dxa0[iobs])
        obs_list.append((t, iobs, obs_val, obs_err_stdv))
    return obs_list
    
# plot: dxa, obs-H(xb), fakedobs-H(xb), for obs_list/faked_obs_list, each item [it,ix,val,err] or [ix,val,err]
def plot_fakedobs_test(xb0, dxa0, obs_list, faked_obs_list):
    nx = xb0.shape[0]
    fig=plt.figure(figsize=(10,6))
    xcoord = np.arange(0,nx)
    zeros = np.zeros((nx,))
    plt.plot(xcoord, zeros, 'b--')
    plt.plot(xcoord, dxa0, 'r', label='xa - xb')
    #get obs-H(xb)
    if(len(obs_list[0]) == 3):
        obs_xs, obs_vals, obs_errs = zip(*obs_list)
    elif(len(obs_list[0]) == 4):
        tmp, obs_xs, obs_vals, obs_errs = zip(*obs_list)
    obs_mb = np.array(obs_vals) - xb0[np.array(obs_xs)]
    #get fakedObs-H(xb)
    if(len(faked_obs_list[0]) == 3):
        obs_xs, obs_vals, obs_errs = zip(*faked_obs_list)
    elif(len(faked_obs_list[0]) == 4):
        tmp, faked_obs_xs, faked_obs_vals, faked_obs_errs = zip(*faked_obs_list)
    faked_obs_mb = np.array(faked_obs_vals) - xb0[np.array(faked_obs_xs)]
    plt.plot(obs_xs, obs_mb, 'g*', label='obs-H(xt)')
    plt.plot(faked_obs_xs, faked_obs_mb, 'c.', label='fakedobs-xt')
    plt.legend()
    plt.xlabel('x coord')
    plt.savefig("fakedobs_test_plot.png")
    plt.close()



if(len(sys.argv[1:]) == 1):
    fileconf = sys.argv[1]
else:
    help()

with open(fileconf, "r") as f:
    dct = json.load(f)["fake_obs"]
    nx=dct["nx"]
    it=dct["it"]
    obs_dict = read_obs_dict(dct["obs_csv"]) 
    xas = read_xs_bin(dct["3dvar_xas_file"], nx)
    xbs = read_xs_bin(dct["xbs_file"], nx)
    dxa0 = xas[it,:] - xbs[it,:]
    
    faked_obs_list = makeFakeObs(dxa0, xbs[it,:], it, dct["slam_xb0_err_stdv"], dct["slam_obs_err_stdv"])
    output_csv(faked_obs_list, dct["output_fakeobs_csv"])
    
    #plot test 
    plot_fakedobs_test(xbs[it,:], dxa0, obs_dict[it], faked_obs_list)

    
