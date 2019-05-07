#!/usr/bin/env python

import os
import numpy as np
import glob

# concat run_dir/*/output.bin's istep to concat_file_path
def concat_ensemble_astep(run_dir, istep, output_file_name, concat_file_path, nx):
    files=glob.glob(run_dir+"ens/ens_*/"+output_file_name)
    files.sort()
    nfile = len(files)
    concat_data = np.zeros((nfile, nx), dtype="double")
    for i, filename in enumerate(files):
        Xs=np.fromfile(filename, dtype="double", count=-1)
        nt = int(len(Xs) / nx)
        Xs=Xs.reshape((nt, nx))
        concat_data[i,:] = Xs[istep,:]
    concat_data.tofile(concat_file_path)

def concat_ensemble_steps(run_dir, steps, output_file_name, concat_file_prefix, nx):
    files=glob.glob(run_dir+"/ens_*/"+output_file_name)
    files.sort()
    nfile = len(files)
    nstep = len(steps)
    concat_data = np.zeros((nstep, nfile, nx), dtype="double")
    for i, filename in enumerate(files):
        Xs=np.fromfile(filename, dtype="double", count=-1)
        nt = int(len(Xs) / nx)
        Xs=Xs.reshape((nt, nx))
        for j, jstep in enumerate(steps):
            concat_data[j, i,:] = Xs[jstep,:]
    for j, jstep in enumerate(steps):
        print(jstep)
        concat_file_path = concat_file_prefix+"%0.5d.bin"%jstep
        concat_data[j,:,:].tofile(concat_file_path)



steps=288
n_ens=10000
cwd=os.getcwd()
template_burgers_dir=cwd+"/burgers_template_288steps/"
run_dir=cwd+"/Q288_fullwindow_10K/"
x_init_file=cwd+"/x_init/x_init_000.bin"
concat_file_path=run_dir+"concat_ens_step288.bin"
'''
for iens in range(n_ens):
    mem_dir=run_dir+"/ens_%0.5d" %iens
    if(os.path.isdir(mem_dir)):
        continue
    os.mkdir(mem_dir)
    os.chdir(mem_dir)
    os.symlink(template_burgers_dir+"/burgers_conf.in", "burgers_conf.in")
    os.symlink(template_burgers_dir+"/run_burgers.x", "run_burgers.x")
    os.symlink(x_init_file,"x_init.bin")
    os.system("./run_burgers.x")
'''
#concat
#concat_ensemble_astep(run_dir, steps, "output.bin", concat_file_path, 1000)
#concat_ensemble_steps(run_dir, np.arange(288,289,1), "output.bin", "Q288_fullwindow_10K/concat_ens_step", 1000)
#concat_ensemble_astep(run_dir, 7, "output.bin", run_dir+"concat_ens_step007.bin", 1000)
#concat_ensemble_astep(run_dir, 8, "output.bin", run_dir+"concat_ens_step008.bin", 1000)
#concat_ensemble_astep(run_dir, 9, "output.bin", run_dir+"concat_ens_step009.bin", 1000)
#concat_ensemble_astep(run_dir, 10, "output.bin", run_dir+"concat_ens_step010.bin", 1000)
#concat_ensemble_astep(run_dir, 11, "output.bin", run_dir+"concat_ens_step011.bin", 1000)
#concat_ensemble_astep(run_dir, 12, "output.bin", run_dir+"concat_ens_step012.bin", 1000)

concat_ensemble_astep(run_dir, 36, "output.bin", run_dir+"concat_ens_step036.bin", 1000)
concat_ensemble_astep(run_dir, 72, "output.bin", run_dir+"concat_ens_step072.bin", 1000)
concat_ensemble_astep(run_dir, 108, "output.bin", run_dir+"concat_ens_step108.bin", 1000)
concat_ensemble_astep(run_dir, 144, "output.bin", run_dir+"concat_ens_step144.bin", 1000)
concat_ensemble_astep(run_dir, 216, "output.bin", run_dir+"concat_ens_step216.bin", 1000)
concat_ensemble_astep(run_dir, 288, "output.bin", run_dir+"concat_ens_step288.bin", 1000)
