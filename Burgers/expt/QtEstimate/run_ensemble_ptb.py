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




steps=6
n_ens=10000
cwd=os.getcwd()
template_burgers_dir=cwd+"/burgers_template"
run_dir=cwd+"/Q006_10K/"
x_init_file=cwd+"/x_init/x_init_000.bin"
concat_file_path=run_dir+"concat_ens_step006.bin"
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
#concat_ensemble_astep(run_dir, 6, "output.bin", concat_file_path, 1000)
