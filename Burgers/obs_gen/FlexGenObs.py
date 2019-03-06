#!/usr/bin/env python
import numpy as np
import pprint as pp
import csv
import json
import sys

#return list: [ (t, x, obs, err) ] 
def add_noise(Xs, obs_ts, obs_xs, obs_error_mean, obs_error_stdv):
    res=[]
    noise = np.random.normal(obs_error_mean, obs_error_stdv, (len(obs_ts), len(obs_xs)))
    it=0
    for t in obs_ts:
        ix=0
        for x in obs_xs:
            res.append((t, x, Xs[t,x] + noise[it,ix], obs_error_stdv))
            ix+=1
        it+=1
    return res
        
def output_csv(res, outfile):
    with open(outfile, "w") as f:
        writer = csv.writer(f)
        writer.writerow(["t","x", "truth", "obs", "err_stdv"])
        for tup in res:
            writer.writerow(tup)
    

def help():
    print( '''
usage: FlexGenObs.py [model state.bin] [fgobs_conf.in] [output obs.csv] [nx, default 1000]
''')
    exit(0)


if(len(sys.argv[1:]) in [3,4] ):
    filebin = sys.argv[1]
    fileconf = sys.argv[2]
    fileout = sys.argv[3]
    if(len(sys.argv[1:]) == 4):
        nx = int(sys.argv[4])
    else:
        nx = 1000
else:
    help()

Xs1=np.fromfile(filebin, dtype="double", count=-1)
Nt = int(len(Xs1) / nx)
Xs1=Xs1.reshape((Nt, nx))

res_obs=[]
with open(fileconf, "r") as f:
    dct = json.load(f)
    cnt=0
    for obs_set in dct["obs_config_array"]:
        cnt+=1
        print("obs set %d --------------" %(cnt))
        print(type(obs_set))
        if(obs_set["obs_x_end"] == -1):
            obs_set["obs_x_end"] = nx-1
        if(obs_set["obs_t_end"] == -1):
            obs_set["obs_t_end"] = Nt-1
        obs_xs = list(range(obs_set["obs_x_start"], obs_set["obs_x_end"]+1, obs_set["obs_x_step"]))
        obs_ts = list(range(obs_set["obs_t_start"], obs_set["obs_t_end"]+1, obs_set["obs_t_step"]))
        print("nx_obs=%d, nt_obs=%d, n_obs=%d" %(len(obs_xs), len(obs_ts), len(obs_xs)*len(obs_ts)))
        res_obs.extend(add_noise(Xs1, obs_ts, obs_xs, obs_set["obs_error_mean"], obs_set["obs_error_stdv"]))


    output_csv(res_obs, fileout)
