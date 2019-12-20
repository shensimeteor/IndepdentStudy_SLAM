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
QuickPlots.py <truth.bin>. <nx, default 1000> (Note, run in slam or 4dvar running dir)
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

#RMSE time-series (inside window) of background (blue), analysis (red)
def plot_dxb_dxa_timeseries(Xs_truth, Xs_bckg, Xs_anal, opt_metric = "rmse"): #{{{
    if(opt_metric == "rmse"):
        bckg_metric=np.sqrt(np.mean((Xs_bckg - Xs_truth)**2, axis=1))
        anal_metric=np.sqrt(np.mean((Xs_anal - Xs_truth)**2, axis=1))
    elif(opt_metric == "stdv"):
        bckg_metric=np.std(Xs_bckg - Xs_truth, axis=1)
        anal_metric=np.std(Xs_anal - Xs_truth, axis=1)
    elif(opt_metric == "bias"):
        bckg_metric=np.mean(Xs_bckg - Xs_truth, axis=1)
        anal_metric=np.mean(Xs_anal - Xs_truth, axis=1)
    else:
        print("Error, opt_metric unknown")
        return
    fig=plt.figure(figsize=(10,4))
    nt,nx=Xs_truth.shape
    tcoord = np.arange(0,nt)
    plt.plot(tcoord, bckg_metric, 'b', label="xb")
    plt.plot(tcoord, anal_metric, 'r', label="xa")
    plt.xlabel("time step")
    plt.ylabel(opt_metric)
    plt.legend()
    plt.savefig("images/%s_timeseries.png" %opt_metric)
    plt.close()
    return (bckg_metric, anal_metric)
#}}}

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
def plot_atime_xbxaobs_err(xts, xbs, xas, obs_dict, t_step=1, xbs_initguess=None):
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
        if (xbs_initguess is not None):
            plt.plot(xcoord, xbs_initguess[t,:] - xts[t,:], 'c', label='xbIG-xt')
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
    if(n_arg == 1):
        xts_file = sys.argv[1]
        nx=1000
    elif(n_arg == 2):
        xts_file = sys.argv[1]
        nx = int(sys.argv[2])
    else:
        help()

    check_file_exist(["obs.csv", "xbs.bin", "xas.bin"])
    
    print("start reading files")
    xbs=np.fromfile("xbs.bin", dtype="double", count=-1)
    n = xbs.shape[0]
    nt = int(n/nx)
    xbs=xbs.reshape((nt,nx))
    
    xas=np.fromfile("xas.bin", dtype="double", count=-1)
    xas=xas.reshape((nt,nx))

    xts=np.fromfile(xts_file, dtype="double", count=-1)
    xts=xts.reshape((nt,nx))

    # xb init
    xbs_ig=np.fromfile("xbs_initguess.bin", dtype="double", count=-1)
    xbs_ig=xbs_ig.reshape((nt,nx))

    obs_dict = read_obs_dict("obs.csv")

    if(not os.path.isdir("images")):
        os.mkdir("images")

    print("start plotting RMSE time-series")
    plot_dxb_dxa_timeseries(xts, xbs, xas)

    print("start plotting Xb/Xa/Xt/Obs")
    plot_atime_xbxaxt_obs(xts, xbs, xas, obs_dict, 36) #36=6*6, every 6 hour one plot
    plot_atime_xbxaobs_err(xts, xbs, xas, obs_dict, 36, xbs_initguess=xbs_ig)

    

