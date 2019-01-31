#!/usr/bin/env python
import numpy as np
import pprint as pp
import csv

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
        writer.writerow(["t","x", "obs", "err_stdv"])
        for tup in res:
            writer.writerow(tup)
    


nx=1000
N=1441  #per 600s
file_burgers="burgers_1e-6_noise.bin"
Xs=np.fromfile(file_burgers, dtype="double", count=-1).reshape(N,nx)
print(Xs.shape)

window_nt=6*48  #48 hour
obs_ts = list(range(0,window_nt+2,12))  #obs every 2hour
print(obs_ts)

obs_xs = list(range(0,nx, 100))
print(obs_xs)

obs_error_mean=0.0
obs_error_stdv=0.5

res=add_noise(Xs, obs_ts, obs_xs, obs_error_mean, obs_error_stdv)
output_csv(res, "burgers_1e-6_noise_obs_sample.csv")

