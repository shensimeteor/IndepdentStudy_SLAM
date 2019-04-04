#!/usr/bin/env python

#run batch of 3dvars (use same B) & fake obs

import os
import csv


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

def obs_3col_to_4ccol(obs_list, t):
    new_list = []
    for obs in obs_list:
        new_list.append([t, obs[0], obs[1], obs[2]])
    return new_list

# output obs csv
def output_obs_csv(res, outfile, t=None):
    with open(outfile, "w") as f:
        writer = csv.writer(f)
        writer.writerow(["t","x", "obs", "err_stdv"])
        for tup in res:
            if(len(tup) == 3 and t is not None):
                item = [t, tup[0], tup[1], tup[2]]
            elif(len(tup) == 4):
                item = tup
            writer.writerow(item)



# return list of faked_obs [(t,x,val,err)]
def prepare_dir_3dvar_fakeobs(srcdir, tgtdir, step, obs_list):
    os.symlink(srcdir+"xbs.bin", "xbs.bin")
    print("in dir %s -----------------" %tgtdir)
    os.system("ExtractAPeriod.py xbs.bin xb0.bin %d 1" %step)
    output_obs_csv(obs_list, "ob0.csv", 0)
    os.symlink(srcdir+"burgers_conf.in", "burgers_conf.in")
    os.symlink(srcdir+"run_burgers_4dvar.x", "run_burgers_4dvar.x")
    os.symlink(srcdir+"cov", "cov")
    #run
    os.system("./run_burgers_4dvar.x > log.run_burgers_4dvar")
    ##run plot
    #os.system("QuickPlotATimeStep.py ../truth.bin")
    #run fakeobs
    os.symlink(srcdir+"fakeobs_conf.json", "fakeobs_conf.json")
    os.system("FakeObsFrom3DVar.py fakeobs_conf.json")


template_3dvar_dir = "/home/si/Work/IndepdentStudy_SLAM//Burgers/expt/setting2/slam_tunes/FakedObs/alltimes_FakeObs/template_3dvar/"
obs_file = "/home/si/Work/IndepdentStudy_SLAM//Burgers/expt/setting2/slam_tunes/FakedObs/alltimes_FakeObs/obs.csv"
xts_file = "/home/si/Work/IndepdentStudy_SLAM//Burgers/expt/setting2/slam_tunes/FakedObs/alltimes_FakeObs/truth.bin"

obs_dict = read_obs_dict(obs_file)

all_faked_obs_list = []

for istep in range(0,289,6):
    workdir="3dvar_fakeobs_%0.3d" %istep
#    os.mkdir(workdir)
    os.chdir(workdir)
#    prepare_dir_3dvar_fakeobs(template_3dvar_dir, workdir, istep, obs_dict[istep])
    faked_obs_dict = read_obs_dict("faked_ob0.csv")
    all_faked_obs_list.extend( obs_3col_to_4ccol(faked_obs_dict[0], istep))
    os.chdir("..")

print(len(all_faked_obs_list))
print(all_faked_obs_list[0])
output_obs_csv(all_faked_obs_list, "faked_obs.csv")
